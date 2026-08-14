# ==========================================
# 0. RUST CHECKS (Runs in Terminal via R)
# ==========================================

# preliminaries to avoid linker error
#####################################
# Define the base path to the GCC libraries in Rtools45
rtools_gcc_base <- "C:/rtools45/x86_64-w64-mingw32.static.posix/lib/gcc/x86_64-w64-mingw32.static.posix"

# Find the versioned folder inside it (e.g., "13.2.0")
gcc_version_folder <- list.files(rtools_gcc_base)[1]
full_lib_path <- file.path(rtools_gcc_base, gcc_version_folder)

# Define file paths
source_file <- file.path(full_lib_path, "libgcc.a")
target_file <- file.path(full_lib_path, "libgcc_eh.a")

# Copy and rename the file
if (file.exists(source_file) && !file.exists(target_file)) {
  file.copy(source_file, target_file)
  message("Success! libgcc_eh.a created.")
} else if (file.exists(target_file)) {
  message("libgcc_eh.a already exists. You are good to go!")
} else {
  warning("Could not find libgcc.a. Please check your Rtools path.")
}
#################################################################

# Run format, clippy, and tests using the GNU toolchain
# Format code
system2("cargo", args = c("fmt", "--manifest-path", "src/rust/Cargo.toml"))

# Run clippy (linter)
system2("cargo", args = c("clippy", "--manifest-path", "src/rust/Cargo.toml"))

# Run Rust unit tests with the GNU toolchain
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
