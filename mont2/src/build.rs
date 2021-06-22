extern crate cc;

fn main() {
    cc::Build::new().file("src/mont_mul_asm.S").compile("mont_mul_asm");
}
