//!
//! Nadgrids multi threaded catalog
//!
//! Maintain a list of loaded grids
//!
use super::grid::Nadgrid;
use std::ptr::null_mut;
use std::sync::atomic::{AtomicPtr, Ordering};
use std::sync::Mutex;

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
    next: AtomicPtr<Node>,
}

impl Node {
    fn new(grid: Nadgrid) -> Self {
        Self {
            grid,
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
}

/// Private catalog implementation
pub(super) struct Catalog {
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
    fn iter(&self) -> impl Iterator<Item = &'static Node> {
        std::iter::successors(Node::get(&self.first), |prev| Node::get(&prev.next))
    }

    /// Add an externally created grid
    /// to the catalog
    fn add_node(&self, grid: Nadgrid) -> &'static Node {
        let last = self.iter().last().map(|n| &n.next).unwrap_or(&self.first);
        let node_ptr = Box::into_raw(Box::new(Node::new(grid)));
        last.store(node_ptr, Ordering::Relaxed);
        unsafe { &*node_ptr }
    }

    /// Find a grid from its name
    fn find(&self, name: &str) -> Option<GridRef> {
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
    fn set_builder(&mut self, builder: GridBuilder) -> Option<GridBuilder> {
        self.builder.replace(builder)
    }

    fn add_grid(&self, grid: Nadgrid) {
        self.add_node(grid);
    }
}

pub(crate) mod catalog {
    use super::*;
    use lazy_static::lazy_static;

    lazy_static! {
        static ref CATALOG: Mutex<Catalog> = Mutex::new(Catalog::default());
    }

    pub(crate) fn find_grid(name: &str) -> Option<GridRef> {
        CATALOG.lock().unwrap().find(name)
    }

    pub(crate) fn add_grid(grid: Nadgrid) {
        CATALOG.lock().unwrap().add_grid(grid)
    }

    pub(crate) fn set_builder(builder: GridBuilder) -> Option<GridBuilder> {
        CATALOG.lock().unwrap().set_builder(builder)
    }
}
