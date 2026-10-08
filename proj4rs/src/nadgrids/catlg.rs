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

#[cfg(not(target_arch = "wasm32"))]
mod implem {
    use super::Node;
    use std::ptr::null_mut;
    use std::sync::atomic::{AtomicPtr, Ordering};

    #[derive(Debug)]
    pub(super) struct NodePtr(AtomicPtr<Node>);

    impl Default for NodePtr {
        #[inline]
        fn default() -> Self {
            Self::new()
        }
    }

    impl NodePtr {
        #[inline]
        pub(super) const fn new() -> Self {
            Self(AtomicPtr::new(null_mut::<Node>()))
        }
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

#[cfg(target_arch = "wasm32")]
mod implem {
    use super::Node;
    use std::cell::Cell;

    #[derive(Debug)]
    pub(super) struct NodePtr(Cell<Option<&'static Node>>);

    impl Default for NodePtr {
        #[inline]
        fn default() -> Self {
            Self::new()
        }
    }

    impl NodePtr {
        #[inline]
        pub(super) const fn new() -> Self {
            Self(Cell::new(None))
        }
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
pub struct Node {
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
}

#[cfg(target_arch = "wasm32")]
mod arch {
    use super::GridBuilder;
    use std::cell::RefCell;
    pub type BuilderRef = RefCell<Option<GridBuilder>>;
    pub const EMPTY_BUILDER_REF: BuilderRef = RefCell::new(None);
}

#[cfg(not(target_arch = "wasm32"))]
mod arch {
    use super::GridBuilder;
    pub type BuilderRef = Option<GridBuilder>;
    pub const EMPTY_BUILDER_REF: BuilderRef = None;
}

pub struct Catalog {
    first: NodePtr,
    builder: arch::BuilderRef,
}

impl Catalog {
    const fn new() -> Self {
        Self {
            first: NodePtr::new(),
            builder: arch::EMPTY_BUILDER_REF,
        }
    }

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
        // All nodes with the same name are from the same file
        let mut iter = self.iter().filter(move |n| n.name == name).peekable();
        iter.peek().is_some().then(|| iter.map(|n| &n.grid))
    }

    /// Add a grid to the gridlist
    /// Note that parent must exists in the list.
    pub fn add_grid(&self, name: String, grid: Grid) -> Result<(), Error> {
        let parent = if !grid.is_root() {
            self.iter()
                .find(|n| n.name == name && n.grid.id == grid.lineage)
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

impl Default for Catalog {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(not(target_arch = "wasm32"))]
pub mod catalog {
    use super::*;
    use std::sync::Mutex;

    static CATALOG: Mutex<Catalog> = Mutex::new(Catalog::new());

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
                    #[allow(unused)]
                    if let Err(err) = b(&cat, name) {
                        error!("Error looking for grid shift {name}: {err}");
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

    pub fn with_catalog<F, R>(f: F) -> R
    where
        F: FnOnce(&Catalog) -> R,
    {
        let ctlg = CATALOG.lock().unwrap();
        f(&ctlg)
    }
}

#[cfg(target_arch = "wasm32")]
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

    pub fn with_catalog<F, R>(f: F) -> R
    where
        F: FnOnce(&Catalog) -> R,
    {
        CATALOG.with(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nadgrids;
    use crate::tests::setup;

    #[test]
    #[cfg(feature = "local_tests")]
    fn test_nadgrid_same_subgrid_names_in_files() {
        // Regression test: the parent of a subgrid must be searched only among the
        // grids of the same file.
        //
        // ntv2_0.gsb and MAY76V20.gsb both have root grids named CAeast and
        // CAwest with subgrids: the subgrids of one file must not be attached
        // to the roots of the other one.
        setup();

        catalog::set_builder(nadgrids::files::read_from_file);

        // Load both files
        for file in ["north-america/ntv2_0.gsb", "north-america/MAY76V20.gsb"] {
            assert!(catalog::find_grids(file, &mut vec![]), "{file} not found");
        }

        // Collect bad links and check outside the catalog lock, so that a failure
        // does not poison the catalog for the other tests.
        let foreign_parents: Vec<_> = catalog::with_catalog(|cat| {
            cat.iter()
                .filter_map(|node| node.parent.map(|parent| (node, parent)))
                .filter(|(node, parent)| parent.name != node.name)
                .map(|(node, parent)| format!("{} -> {}", node.name, parent.name))
                .collect()
        });

        assert!(
            foreign_parents.is_empty(),
            "subgrids attached to a parent from another file: {foreign_parents:?}",
        );
    }
}
