use std::env;
use std::path::Path;

fn main() {
    // 1. Tell cargo to link the 'RNA' library (libRNA.so / libRNA.dylib)
    println!("cargo:rustc-link-lib=RNA");

    // 2. Allow users to override the search path via an environment variable
    if let Ok(custom_dir) = env::var("VIENNARNA_LIB_DIR") {
        println!("cargo:rustc-link-search=native={}", custom_dir);
        return;
    }

    // 3. List standard default library search paths across different OS architectures
    let default_paths = vec![
        "/opt/homebrew/lib", // macOS Apple Silicon
        "/usr/local/lib",    // macOS Intel & general Linux builds from source
        "/usr/lib",          // Standard Linux package managers
        "/usr/lib/x86_64-linux-gnu", // Ubuntu/Debian 64-bit
    ];

    for path in default_paths {
        if Path::new(path).exists() {
            println!("cargo:rustc-link-search=native={}", path);
        }
    }
}