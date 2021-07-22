# Comments on initial CRAN submission

July 19, 2021

## Test environments

- R version 4.0.5 (2021-03-31), Platform: x86_64-w64-mingw32/x64 (64-bit), Running under: Windows 10 x64 (build 17134)
- R version 4.0.3 (2020-10-10), Platform: x86_64-pc-linux-gnu (64-bit), Running under: Ubuntu 16.04.7 LTS

## R CMD check results

No ERRORs, WARNINGs, or NOTEs on either platform.

# Comments on second CRAN submission

July 20, 2021

This submission corrects issues flagged in the initial submission of July 19.

- License changed to MIT and LICENSE file formatted accordingly.
- URLs in documentation corrected.
- Title case for package title.

## Test environments

- R version 4.0.5 (2021-03-31), Platform: x86_64-w64-mingw32/x64 (64-bit), Running under: Windows 10 x64 (build 17134)
- R version 4.0.3 (2020-10-10), Platform: x86_64-pc-linux-gnu (64-bit), Running under: Ubuntu 16.04.7 LTS

## R CMD check results

No ERRORs, WARNINGs, or NOTEs on either platform.

# Comments on third CRAN submission

July 22, 2021

This submission corrects issues raised by Uwe and Gregor from the previous submission.

- Title case corrected.
- DESCRIPTION no longer starts with "The Ostats package ..."
- Example is no longer wrapped in dontrun.
- Argument verbose added to Ostats() and Ostats_multivariate(). If set to TRUE, message() is used to display progress messages, rather than print().
- References added to the description field of the DESCRIPTION file. I only included references to our group's own published work in the DESCRIPTION file. Other references listed in the README.md file are for further reading and not appropriate for the DESCRIPTION file.
- Added BugReports field to DESCRIPTION.
