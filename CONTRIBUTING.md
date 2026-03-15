# Contributing

Thank you for improving this research software project.

## Contribution Scope

Contributions are welcome for:

- simulation correctness and performance
- PBPK/ABM methodology improvements
- analysis and reproducibility tooling
- documentation quality
- test coverage and regression prevention

## Development Workflow

1. Fork and create a feature branch.
2. Keep commits focused and descriptive.
3. Add or update tests when behavior changes.
4. Update docs for new scripts, arguments, or outputs.
5. Open a pull request with a concise scientific and technical summary.

## Pull Request Checklist

- [ ] The change is reproducible and has a clear motivation.
- [ ] Commands, scripts, and file paths in docs are correct.
- [ ] New outputs or metrics are documented.
- [ ] Any assumptions, limitations, and risks are stated.
- [ ] Relevant citations are preserved or updated.

## Coding Guidance

- Prefer deterministic defaults where possible.
- Keep performance-sensitive loops explicit and benchmarkable.
- Avoid hard-coded machine-specific paths.
- Preserve compatibility with existing output parsing workflows.

## Issue Reporting

Please include:

- operating system and tool versions (compiler, MATLAB, Python)
- exact command(s) used
- input config or dosing file
- expected behavior vs observed behavior
- relevant log lines or error trace

## Scientific Integrity

When proposing model changes:

- explain biological rationale
- identify affected outputs/figures
- state whether previous results are expected to shift
- provide a minimal reproduction case

## License

By contributing, you agree that your contributions will be licensed under the BSD 3-Clause License in this repository.
