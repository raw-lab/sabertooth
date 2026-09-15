//! ASCII art banners.
//!
//! One sabertooth mascot, four framings — shown by `version`, `info`, `doctor`
//! and `help`. Kept as plain ASCII so it renders in any terminal / log capture.

/// The core sabertooth skull-and-fangs mascot.
pub const SABERTOOTH: &str = r#"
                     /\   /\
                    //\\_//\\     ____
                    \_     _/    /   /
                     / * * \    /^^^]
                     \_\O/_/    [   ]
                      /   \_    [   /
                      \     \_  /  /
                       [ [ /  \/ _/
                      _[ [ \  /_/
                     ( \_\_\/ /
                      \  V V /   S A B E R T O O T H
                       \_\_\/
                      /  /\  \
                     /_/  \_\_\
                    ((    ))       ,,,,,
                     \\  //       /     \
                    __\\//__     | fangs |
                   /   ||   \     \_____/
"#;

/// Compact fang-mark used inline in headers.
pub const FANGS: &str = r#"   \V/   \V/ "#;

/// Banner shown by `version`.
pub fn version_banner(version: &str) -> String {
    format!(
        r#"
   ___       _           _              _   _
  / __| __ _| |__  ___ _| |_ ___  ___ _| |_| |_
  \__ \/ _` | '_ \/ -_)  _/ _ \/ _ \  _|  _| ' \
  |___/\__,_|_.__/\___|\__\___/\___/\__|\__|_||_|
       \V/   \V/    fast · light · profile-sensitive

  Sabertooth v{version}
  A pure-Rust reimplementation of the MMseqs2 profile/PSSM core.
"#,
        version = version
    )
}

/// Banner shown by `info`.
pub const INFO_BANNER: &str = r#"
   ┌───────────────────────────────────────────────┐
   │   \V/        S A B E R T O O T H        \V/    │
   │    ▔▔    profile-sensitive sequence search     │
   └───────────────────────────────────────────────┘
"#;

/// Banner shown by `doctor`.
pub const DOCTOR_BANNER: &str = r#"
      \V/   \V/
   ___ \\___// ___     S A B E R T O O T H  ::  doctor
  /   \ (o o) /   \    self-check & environment report
  \___/  \_/  \___/
       fangs sharp?
"#;

/// Banner shown by `help`.
pub const HELP_BANNER: &str = r#"
   \V/  S A B E R T O O T H  \V/
  ══════════════════════════════
   fast · light · profile-sensitive
"#;
