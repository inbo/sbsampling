#![recursion_limit = "256"]
use extendr_api::prelude::*;

mod rng;
mod sample_hex;
mod sobol;

// ---------------------------------------------------------------------------
// extendr module registration
// ---------------------------------------------------------------------------

extendr_module! {
    mod sbsampling;
    use sample_hex;
    use sobol;
}
