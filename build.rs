extern crate git_build_version;

const PACKAGE_TOP_DIR: &'static str = ".";

fn main() {
    git_build_version::write_version(PACKAGE_TOP_DIR).expect("Saving git version");
}
