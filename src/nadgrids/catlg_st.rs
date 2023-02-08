//!
//! Nadgrids single threaded catalog
//!
//! Maintain a list of loaded grids
//!
use super::grid::Grid;
use crate::errors::Error;
use std::cell::{Cell, RefCell};

/// Nadgrid factory: simple function pointer that return a NadGrid.
///
/// This is an infaillible method that should return [`None`] if
/// no Nadgrid can be found or if an error occured when loading/building
/// the nadgrid.
pub(crate) type GridBuilder = fn(&Catalog, &str) -> Result<(), Error>;

/// Static reference to nadgrids
///
/// Grids  have a static lifetime on the heap
/// It means they are never deallocated;
pub(crate) type GridRef = &'static Grid;

/// Node to chain loaded nadgrids
struct Node {
    name: String,
    grid: Grid,
    parent: Option<&'static Node>,
    next: Cell<Option<&'static Node>>,
}

impl Node {
    fn new(name: String, grid: Grid, parent: Option<&'static Node>) -> Self {
        Self {
            name,
            grid,
            parent,
            next: Cell::new(None),
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
#[derive(Default)]
pub(crate) struct Catalog {
    first: Cell<Option<&'static Node>>,
    builder: RefCell<Option<GridBuilder>>,
}

impl Catalog {
    /// Set a builder callback, None if no builder
    /// was set.
    fn set_builder(&self, builder: GridBuilder) -> Option<GridBuilder> {
        self.builder.borrow_mut().replace(builder)
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
        let node = Box::leak::<'static>(Box::new(node));
        if let Some(parent) = node.parent {
            // Insert the node juste behind parent
            node.next.replace(parent.next.replace(Some(node)));
        } else {
            let last = self.iter().last().map(|n| &n.next).unwrap_or(&self.first);
            last.replace(Some(node));
        }
        node
    }

    pub(crate) fn find(&self, name: &str) -> Option<impl Iterator<Item = GridRef>> {
        let mut iter = self.iter();
        let node = iter.find(|n| n.name == name);
        node.map(|node| {
            std::iter::once(&node.grid).chain(iter.filter(|n| n.is_child_of(node)).map(|n| &n.grid))
        })
    }

    /// Add a grid to the gridlist
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

    thread_local! {
        static CATALOG: Catalog = Catalog::default();
    }

    pub(crate) fn find_grids(name: &str, grids: &mut Vec<GridRef>) -> bool {
        CATALOG.with(|cat| match cat.find(name) {
            Some(iter) => {
                grids.extend(iter);
                true
            }
            None => cat
                .builder
                .borrow()
                .and_then(|b| {
                    b(cat, name);
                    cat.find(name).map(|iter| grids.extend(iter))
                })
                .is_some(),
        })
    }

    pub(crate) fn add_grid(name: String, grid: Grid) -> Result<(), Error> {
        CATALOG.with(|cat| cat.add_grid(name, grid))
    }

    pub(crate) fn set_builder(builder: GridBuilder) -> Option<GridBuilder> {
        CATALOG.with(|cat| cat.set_builder(builder))
    }
}
