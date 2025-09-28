[![Build + Test](https://github.com/LBNL-ETA/Radiance/actions/workflows/build.yaml/badge.svg)](https://github.com/LBNL-ETA/Radiance/actions/workflows/build.yaml)
[![GitHub All Releases](https://img.shields.io/github/downloads/LBNL-ETA/Radiance/total?label=Download%20Release)](https://github.com/LBNL-ETA/Radiance/releases)

# Radiance

A Validated Lighting Simulation Tool

This repository is a mirror of the official Radiance CVS source tree from [radiance-online.org](http://www.radiance-online.org). It is automatically updated every day at 00:00 UTC.

## Source

The source code is located on the [`master` branch](https://github.com/LBNL-ETA/Radiance/tree/master).

**Important:** Do not commit directly to the `master` branch. This branch is a one-way sync from the official CVS repository. Any manual edits will be overwritten by the daily import process.

## Testing

Unit tests are run automatically for every build on Windows, macOS (x86-64 and arm64), and Linux. You can view the test suite in the [test directory](https://github.com/LBNL-ETA/Radiance/tree/master/test).

## Installer

Installers for Windows, macOS (x86-64 and arm64), and Linux are built and published with each new commit to `master`.

You can download the latest installers from the [**Releases** page](https://github.com/LBNL-ETA/Radiance/releases). Installers are listed under the "Assets" section of each release.
