#![recursion_limit = "256"]
use extendr_api::prelude::*;

mod grts_h3;
mod sobol;

// ---------------------------------------------------------------------------
// extendr module registration
// ---------------------------------------------------------------------------

extendr_module! {
    mod sbsampling;
    use grts_h3;
    use sobol;
}
