//! Process measurements (peak RSS, page faults, CPU time) and a tiny JSON
//! builder for machine-readable reports (no serde in the default build).

/// Snapshot of process resource usage (Linux `/proc`; zeros elsewhere).
#[derive(Clone, Copy, Debug, Default)]
pub struct ProcStats {
    /// Peak resident set size (`VmHWM`), bytes.
    pub peak_rss: u64,
    /// Current resident set size (`VmRSS`), bytes.
    pub rss: u64,
    /// Minor / major page faults.
    pub minflt: u64,
    pub majflt: u64,
    /// User / system CPU seconds.
    pub utime: f64,
    pub stime: f64,
}

/// Read the current process statistics.
pub fn proc_stats() -> ProcStats {
    let mut p = ProcStats::default();
    if let Ok(s) = std::fs::read_to_string("/proc/self/status") {
        for line in s.lines() {
            let kb = |l: &str| -> u64 {
                l.split_whitespace()
                    .nth(1)
                    .and_then(|v| v.parse::<u64>().ok())
                    .unwrap_or(0)
                    * 1024
            };
            if line.starts_with("VmHWM:") {
                p.peak_rss = kb(line);
            } else if line.starts_with("VmRSS:") {
                p.rss = kb(line);
            }
        }
    }
    if let Ok(s) = std::fs::read_to_string("/proc/self/stat") {
        // fields after the parenthesised command name
        if let Some(close) = s.rfind(')') {
            let f: Vec<&str> = s[close + 2..].split_whitespace().collect();
            // f[0] = state (field 3); minflt = field 10, majflt = 12, utime 14, stime 15
            let get = |i: usize| -> u64 { f.get(i).and_then(|v| v.parse().ok()).unwrap_or(0) };
            p.minflt = get(7);
            p.majflt = get(9);
            let hz = 100.0; // USER_HZ on Linux
            p.utime = get(11) as f64 / hz;
            p.stime = get(12) as f64 / hz;
        }
    }
    p
}

/// Minimal JSON object builder.
#[derive(Default)]
pub struct Json {
    parts: Vec<String>,
}

fn esc(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\t' => o.push_str("\\t"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o
}

impl Json {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn num<T: std::fmt::Display>(mut self, k: &str, v: T) -> Self {
        let s = v.to_string();
        let s = if s == "NaN" || s == "inf" || s == "-inf" { "null".to_string() } else { s };
        self.parts.push(format!("\"{}\": {}", esc(k), s));
        self
    }
    pub fn float(self, k: &str, v: f64) -> Self {
        if v.is_finite() {
            self.num(k, format!("{:.6}", v))
        } else {
            self.raw(k, "null")
        }
    }
    pub fn str(mut self, k: &str, v: &str) -> Self {
        self.parts.push(format!("\"{}\": \"{}\"", esc(k), esc(v)));
        self
    }
    pub fn bool(mut self, k: &str, v: bool) -> Self {
        self.parts.push(format!("\"{}\": {}", esc(k), v));
        self
    }
    pub fn raw(mut self, k: &str, v: &str) -> Self {
        self.parts.push(format!("\"{}\": {}", esc(k), v));
        self
    }
    pub fn obj(self, k: &str, v: Json) -> Self {
        let s = v.finish();
        self.raw(k, &s)
    }
    pub fn finish(self) -> String {
        format!("{{{}}}", self.parts.join(", "))
    }
}

/// Parse a byte size such as `1048576`, `512K`, `256M`, `2G`, `1.5G`.
pub fn parse_bytes(s: &str) -> Result<usize, String> {
    let t = s.trim();
    if t.is_empty() {
        return Err("empty size".into());
    }
    let (num, mult) = match t.chars().last().unwrap().to_ascii_uppercase() {
        'K' => (&t[..t.len() - 1], 1u64 << 10),
        'M' => (&t[..t.len() - 1], 1u64 << 20),
        'G' => (&t[..t.len() - 1], 1u64 << 30),
        'T' => (&t[..t.len() - 1], 1u64 << 40),
        'B' => (&t[..t.len() - 1], 1),
        _ => (t, 1),
    };
    let v: f64 = num
        .trim()
        .parse()
        .map_err(|_| format!("invalid size '{}' (use e.g. 512M, 2G)", s))?;
    if v.is_nan() || v <= 0.0 {
        return Err(format!("size must be positive: '{}'", s));
    }
    Ok((v * mult as f64) as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizes_parse() {
        assert_eq!(parse_bytes("1024").unwrap(), 1024);
        assert_eq!(parse_bytes("2K").unwrap(), 2048);
        assert_eq!(parse_bytes("1.5G").unwrap(), 3 << 29);
        assert!(parse_bytes("abc").is_err());
        assert!(parse_bytes("-1G").is_err());
    }

    #[test]
    fn json_builds() {
        let j = Json::new().num("a", 1).str("b", "x\"y").obj("c", Json::new().bool("d", true)).finish();
        assert_eq!(j, "{\"a\": 1, \"b\": \"x\\\"y\", \"c\": {\"d\": true}}");
    }

    #[test]
    fn proc_stats_reads_something_on_linux() {
        let p = proc_stats();
        if cfg!(target_os = "linux") {
            assert!(p.peak_rss > 0 && p.rss > 0);
        }
    }
}
