
# Run format, clippy, and tests using the GNU toolchain
system2("cargo", args = c("fmt", "--manifest-path", "src/rust/Cargo.toml"))

# Run clippy (linter)
rextendr::cargo("clippy")

# Run Rust unit tests with the GNU toolchain

rextendr::cargo("test")

pkgbuild::with_build_tools({
  withr::with_envvar(c("CARGO_TARGET_X86_64_PC_WINDOWS_GNU_LINKER" = "gcc"), {
    system2("cargo",
            args = c("test",
                     "--manifest-path", "src/rust/Cargo.toml",
                     "--target", "x86_64-pc-windows-gnu"))
  })
})

# ==========================================
# 1. R PACKAGE DEVELOPMENT
# ==========================================
# Clean up old build artifacts (crucial after a linker error)
rextendr::clean()

# Regenerate the extendr wrappers and documentation
devtools::document()

# Run R-level unit tests
devtools::test()

# ==========================================
# 2. COMPILE & LOAD
# ==========================================
# Restart the R session to ensure no old DLLs are locked in memory
rstudioapi::restartSession()

# Build and install the package locally (Triggers Cargo --release)
devtools::install()

# (Alternative) devtools::load_all() # runs Rust code in debug mode - slower
