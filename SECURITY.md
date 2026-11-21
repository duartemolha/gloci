# Security policy

## Supported versions

gloci is before 1.0. Security fixes go into the next release; only the latest release is
supported.

| Version | Supported |
|---|---|
| latest release | yes |
| older releases | no: upgrade to the latest |

## Reporting a vulnerability

Please do not open a public issue for a security problem. Report it privately instead:
open the repository's **Security** tab and choose **Report a vulnerability**
(<https://github.com/duartemolha/gloci/security/advisories/new>). Only the maintainers see the
report.

Please include:

- the gloci version (`gloci --version`) and how it was installed (release archive, `cargo install`,
  `pip`, bioconda, from source);
- the command, or the Python code, that shows the problem;
- a small input file that triggers it, if one is needed (attach it to the report, not to a public
  issue);
- what happens, and what you expected.

The maintainers aim to acknowledge a report within a week. A confirmed problem is fixed in a new
release and announced in a GitHub security advisory, crediting the reporter unless they prefer
otherwise.

## What counts

gloci reads files that may come from anywhere, so the most likely problems are inputs that make it
misbehave. Report privately:

- an input that makes gloci read or write memory it should not, or crash in a way that suggests
  memory corruption;
- an input or option that makes gloci write outside the files it was asked to write;
- a small input that makes gloci use an unbounded amount of memory, disk or time, where the tool it
  replaces does not;
- a problem in the release files themselves (archives, wheels, the container image).

An error message, an ordinary crash on a malformed file, or an output that differs from bedtools,
bcftools or samtools is a bug, not a vulnerability: please report it in a public
[issue](https://github.com/duartemolha/gloci/issues/new/choose).

A vulnerability in one of gloci's dependencies is best reported to that project. Dependabot watches
gloci's dependencies, and a fixed version is taken up in the next gloci release.
