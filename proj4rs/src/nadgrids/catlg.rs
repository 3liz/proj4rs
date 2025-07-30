//!
//! Nadgrids single threaded catalog
//!
//! Maintain a list of loaded grids
//!
use super::grid::Grid;
use crate::errors::Error;
use crate::log::error;

/// Nadgrid factory: function pointer that load
/// nadgrid into the catalog
///
/// Should return an error if no Nadgrid can be found or
/// an error occurred when loading/building the nadgrid.
pub type GridBuilder = fn(&Catalog, &str) -> Result<(), Error>;

/// Static reference to nadgrids
///
/// Grids  have a static lifetime on the heap
/// It means they are never deallocated;
#[doc(hidden)]
pub type GridRef = &'static Grid;

#[cfg(feature = "multi-thread")]
mod implem {
    use super::Node;
    use std::ptr::null_mut;
    use std::sync::atomic::{AtomicPtr, Ordering};

    #[derive(Debug)]
    pub(super) struct NodePtr(AtomicPtr<Node>);

    impl Default for NodePtr {
        #[inline]
        fn default() -> Self {
            Self(null_mut::<Node>().into())
        }
    }

    impl NodePtr {
        /// Convert raw ptr to static reference
        pub(super) fn get(&self) -> Option<&'static Node> {
            let p = self.0.load(Ordering::Relaxed);
            unsafe { p.as_ref() }
        }
        pub(super) fn insert(&self, node: Node) -> &'static Node {
            node.next
                .0
                .store(self.0.load(Ordering::Relaxed), Ordering::Relaxed);
            let p = Box::into_raw(Box::new(node));
            self.0.store(p, Ordering::Relaxed);
            unsafe { &*p }
        }
    }
}

#[cfg(not(feature = "multi-thread"))]
mod implem {
    use super::Node;
    use std::cell::Cell;

    #[derive(Debug)]
    pub(super) struct NodePtr(Cell<Option<&'static Node>>);

    impl Default for NodePtr {
        #[inline]
        fn default() -> Self {
            Self(Cell::new(None))
        }
    }

    impl NodePtr {
        /// Convert raw ptr to static reference
        #[inline]
        pub(super) fn get(&self) -> Option<&'static Node> {
            self.0.get()
        }
        pub(super) fn insert(&self, node: Node) -> &'static Node {
            let node = Box::leak::<'static>(Box::new(node));
            node.next.0.replace(self.0.replace(Some(node)));
            node
        }
    }
}

use implem::NodePtr;

/// Node to chain loaded nadgrids
#[derive(Debug)]
struct Node {
    name: String,
    grid: Grid,
    parent: Option<&'static Node>,
    next: NodePtr,
}

impl Node {
    fn new(name: String, grid: Grid, parent: Option<&'static Node>) -> Self {
        Self {
            name,
            grid,
            parent,
            next: NodePtr::default(),
        }
    }

    pub fn is_child_of(&self, node: &Self) -> bool {
        match self.parent {
            Some(p) => std::ptr::eq(p, node) || p.is_child_of(node),
            _ => false,
        }
    }
}

/// Private catalog implementation
#[cfg(not(feature = "multi-thread"))]
use std::cell::RefCell;

#[cfg(not(feature = "multi-thread"))]
type BuilderRef = RefCell<Option<GridBuilder>>;

#[cfg(feature = "multi-thread")]
type BuilderRef = Option<GridBuilder>;

#[derive(Default)]
pub struct Catalog {
    first: NodePtr,
    builder: BuilderRef,
}

impl Catalog {
    fn iter(&self) -> impl Iterator<Item = &'static Node> {
        std::iter::successors(self.first.get(), |prev| prev.next.get())
    }

    /// Add an externally created grid
    /// to the catalog
    ///
    /// The insertion ensure that all child nodes are just behind their
    /// parent node
    fn add_node(&self, node: Node) -> &'static Node {
        (if let Some(parent) = node.parent {
            &parent.next
        } else {
            self.iter().last().map(|n| &n.next).unwrap_or(&self.first)
        })
        .insert(node)
    }

    pub fn find(&self, name: &str) -> Option<impl Iterator<Item = GridRef>> {
        let mut iter = self.iter();
        let node = iter.find(|n| n.name == name);
        node.map(|node| {
            std::iter::once(&node.grid).chain(iter.filter(|n| n.is_child_of(node)).map(|n| &n.grid))
        })
    }

    /// Add a grid to the gridlist
    /// Note that parent must exists in the list.
    pub fn add_grid(&self, name: String, grid: Grid) -> Result<(), Error> {
        let parent = if !grid.is_root() {
            self.iter().find(|n| n.grid.id == grid.lineage)
        } else {
            None
        };
        if !grid.is_root() && parent.is_none() {
            return Err(Error::NadGridParentNotFound);
        }
        self.add_node(Node::new(name, grid, parent));
        Ok(())
    }
}

#[cfg(feature = "multi-thread")]
pub mod catalog {
    use super::*;
    use lazy_static::lazy_static;
    use std::sync::Mutex;

    lazy_static! {
        static ref CATALOG: Mutex<Catalog> = Mutex::new(Catalog::default());
    }

    pub fn find_grids(name: &str, grids: &mut Vec<GridRef>) -> bool {
        let cat = CATALOG.lock().unwrap();
        match cat.find(name) {
            Some(iter) => {
                grids.extend(iter);
                true
            }
            None => cat
                .builder
                .and_then(|b| {
                    if b(&cat, name).is_err() {
                        error!("Error looking for grid shift {}", name);
                    }
                    cat.find(name).map(|iter| grids.extend(iter))
                })
                .is_some(),
        }
    }

    pub fn add_grid(name: String, grid: Grid) -> Result<(), Error> {
        CATALOG.lock().unwrap().add_grid(name, grid)
    }

    pub fn set_builder(builder: GridBuilder) -> Option<GridBuilder> {
        CATALOG.lock().unwrap().builder.replace(builder)
    }
}
#[cfg(not(feature = "multi-thread"))]
pub mod catalog {
    use super::*;

    thread_local! {
        static CATALOG: Catalog = Catalog::default();
    }

    pub fn find_grids(name: &str, grids: &mut Vec<GridRef>) -> bool {
        CATALOG.with(|cat| match cat.find(name) {
            Some(iter) => {
                grids.extend(iter);
                true
            }
            None => cat
                .builder
                .borrow()
                .and_then(|b| {
                    if b(cat, name).is_err() {
                        error!("Error looking for grid shift {}", name);
                    }
                    cat.find(name).map(|iter| grids.extend(iter))
                })
                .is_some(),
        })
    }

    pub fn add_grid(name: String, grid: Grid) -> Result<(), Error> {
        CATALOG.with(|cat| cat.add_grid(name, grid))
    }

    pub fn set_builder(builder: GridBuilder) -> Option<GridBuilder> {
        CATALOG.with(|cat| cat.builder.borrow_mut().replace(builder))
    }
}
