//!
//! Nadgrids multi threaded catalog
//!
//! Maintain a list of loaded grids
//!
use super::grid::Grid;
use crate::errors::Error;
use std::ptr::null_mut;
use std::sync::atomic::{AtomicPtr, Ordering};
use std::sync::Mutex;

/// Nadgrid factory: simple function pointer that return a NadGrid.
///
/// This is an infaillible method that should return [`None`] if
/// no Nadgrid can be found or if an error occured when loading/building
/// the nadgrid.
pub type GridBuilder = fn(&Catalog, &str) -> Result<(), Error>;

/// Static reference to nadgrids
///
/// Nadgrids have a static lifetime on the heap
/// It means they are never deallocated;
pub type GridRef = &'static Grid;

/// Node to chain loaded nadgrids
struct Node {
    name: String,
    grid: Grid,
    parent: Option<&'static Node>,
    next: AtomicPtr<Node>,
}

impl Node {
    fn new(name: String, grid: Grid, parent: Option<&'static Node>) -> Self {
        Self {
            name,
            grid,
            parent,
            next: null_mut::<Node>().into(),
        }
    }

    /// Convert raw ptr to static reference
    fn get(p: &AtomicPtr<Node>) -> Option<&'static Node> {
        let p = p.load(Ordering::Relaxed);
        if p.is_null() {
            None
        } else {
            unsafe { Some(&*p) }
        }
    }

    fn is_child_of(&self, node: &Self) -> bool {
        match self.parent {
            Some(p) => std::ptr::eq(p, node) || p.is_child_of(node),
            _ => false,
        }
    }
}

/// Private catalog implementation
pub(crate) struct Catalog {
    first: AtomicPtr<Node>,
    builder: Option<GridBuilder>,
}

impl Default for Catalog {
    fn default() -> Self {
        Self {
            first: null_mut::<Node>().into(),
            builder: None,
        }
    }
}

impl Catalog {
    /// Set a builder callback, None if no builder
    /// was set.
    fn set_builder(&mut self, builder: GridBuilder) -> Option<GridBuilder> {
        self.builder.replace(builder)
    }

    /// Add an externally created grid
    /// to the catalog
    fn add_node(&self, node: Node) -> &'static Node {
        let node_ptr = if let Some(parent) = node.parent {
            // Insert the node juste behind parent
            node.next
                .store(parent.next.load(Ordering::Relaxed), Ordering::Relaxed);
            let node_ptr = Box::into_raw(Box::new(node));
            parent.next.store(node_ptr, Ordering::Relaxed);
            node_ptr
        } else {
            let node_ptr = Box::into_raw(Box::new(node));
            let last = self.iter().last().map(|n| &n.next).unwrap_or(&self.first);
            last.store(node_ptr, Ordering::Relaxed);
            node_ptr
        };
        unsafe { &*node_ptr }
    }

    fn iter(&self) -> impl Iterator<Item = &'static Node> {
        std::iter::successors(Node::get(&self.first), |prev| Node::get(&prev.next))
    }

    /// Find a grid from its name
    pub(crate) fn find(&self, name: &str) -> Option<impl Iterator<Item = GridRef>> {
        let mut iter = self.iter();
        let node = iter.find(|n| n.name == name);
        node.map(|node| {
            std::iter::once(&node.grid).chain(iter.filter(|n| n.is_child_of(node)).map(|n| &n.grid))
        })
    }

    /// Add a grid to the gridlist
    ///
    /// Note that parent must exists in the list.
    pub(crate) fn add_grid(&self, name: String, grid: Grid) -> Result<(), Error> {
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

pub(crate) mod catalog {
    use super::*;
    use lazy_static::lazy_static;

    lazy_static! {
        static ref CATALOG: Mutex<Catalog> = Mutex::new(Catalog::default());
    }

    pub(crate) fn find_grids(name: &str, grids: &mut Vec<GridRef>) -> bool {
        let cat = CATALOG.lock().unwrap();
        match cat.find(name) {
            Some(iter) => {
                grids.extend(iter);
                true
            }
            None => cat
                .builder
                .and_then(|b| {
                    b(&cat, name);
                    cat.find(name).map(|iter| grids.extend(iter))
                })
                .is_some(),
        }
    }

    pub(crate) fn add_grid(name: String, grid: Grid) -> Result<(), Error> {
        CATALOG.lock().unwrap().add_grid(name, grid)
    }

    pub(crate) fn set_builder(builder: GridBuilder) -> Option<GridBuilder> {
        CATALOG.lock().unwrap().set_builder(builder)
    }
}
