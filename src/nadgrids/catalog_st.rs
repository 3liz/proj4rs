//!
//! Nadgrids single threaded catalog
//!
//! Maintain a list of loaded grids
//!
use super::grid::Nadgrid;
use std::cell::Cell;

/// Nadgrid factory: simple function pointer that return a NadGrid.
///
/// This is an infaillible method that should return [`None`] if
/// no Nadgrid can be found or if an error occured when loading/building
/// the nadgrid.
pub(crate) type GridBuilder = fn(&str) -> Option<Nadgrid>;

/// Static reference to nadgrids
///
/// Nadgrids have a static lifetime on the heap
/// It means they are never deallocated;
pub(crate) type GridRef = &'static Nadgrid;

/// Node to chain loaded nadgrids
struct Node {
    grid: Nadgrid,
    next: Cell<Option<&'static Node>>,
}

impl Node {
    fn new(grid: Nadgrid) -> Self {
        Self {
            grid,
            next: Cell::new(None),
        }
    }
}

/// Private catalog implementation
#[derive(Default)]
pub(super) struct Catalog {
    first: Cell<Option<&'static Node>>,
    builder: Option<GridBuilder>,
}

impl Catalog {
    fn iter(&self) -> impl Iterator<Item = &'static Node> {
        std::iter::successors(self.first.get(), |prev| prev.next.get())
    }

    /// Add an externally created grid
    /// to the catalog
    fn add_node(&self, grid: Nadgrid) -> &'static Node {
        let node = Box::leak::<'static>(Box::new(Node::new(grid)));
        let last = self.iter().last().map(|n| &n.next).unwrap_or(&self.first);
        if last.get().is_none() {
            last.replace(Some(node));
        }
        node
    }

    pub(super) fn find(&self, name: &str) -> Option<GridRef> {
        match self.iter().find(|n| n.grid.name() == name) {
            Some(n) => Some(&n.grid),
            None => self
                .builder
                .and_then(|b| b(name))
                .map(|grid| &self.add_node(grid).grid),
        }
    }

    /// Set a builder callback, None if no builder
    /// was set.
    pub(super) fn set_builder(&mut self, builder: GridBuilder) -> Option<GridBuilder> {
        self.builder.replace(builder)
    }

    pub(super) fn add_grid(&self, grid: Nadgrid) {
        self.add_node(grid);
    }
}

pub(crate) mod catalog {
    use super::*;
    use std::cell::RefCell;

    thread_local! {
        static CATALOG: RefCell<Catalog> = RefCell::new(Catalog::default());
    }

    pub(crate) fn find_grid(name: &str) -> Option<GridRef> {
        CATALOG.with(|cat| cat.borrow().find(name))
    }

    pub(crate) fn add_grid(grid: Nadgrid) {
        CATALOG.with(|cat| cat.borrow().add_grid(grid))
    }

    pub(crate) fn set_builder(builder: GridBuilder) -> Option<GridBuilder> {
        CATALOG.with(|cat| cat.borrow_mut().set_builder(builder))
    }
}
