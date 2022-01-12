# MBSinBLG

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fernandopenaranda.github.io/MBSinBLG.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fernandopenaranda.github.io/MBSinBLG.jl/dev)
[![Build Status](https://github.com/fernandopenaranda/MBSinBLG.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fernandopenaranda/MBSinBLG.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/fernandopenaranda/MBSinBLG.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fernandopenaranda/MBSinBLG.jl)

For the sake of reproducibility, we show in this repo: the code, the data files, and the plotting functions used to generate the figures in our work: ref.

Remarks:  

0. All the code provided here is written in Julia (v1.0 >)
1. Our code uses the syntax of [Quantica.jl](https://github.com/pablosanjose/Quantica)
2. The data files used to build each figure can be either directly accessed at `MBSinBLG/data` or generated as we illustrate in the following:

## Installation

In a Julia terminal: 
```
using Pkg
Pkg.add("https://github.com/fernandopenaranda/MBSinBLG")
```

## Usage

Generate all data corresponding to each figure and store it in a folder named `/data`
```
include("~/MBSinBLG/compute/computefigures.jl")
```
The specifications of the machines used in our calculations are shown in said file.

## Plot

We also provide in `~/MBSinBLG/compute/plotfunctions.jl` the plotting code that generates the figures of our manuscript.
The user just need to specify the path to the .csv file and use the corresponding plotting function, for instance, figure 7 is:
```
using CairoMakie, VegaLite, LaTeXStrings, Colors
using ElectronDisplay
using MBSinBLG
pathfig7 = "MBSinBLG/data/fig7"
fig7plot(pathfig7)

```

## Exported API (additional functions)

In addition, we explicitly export a set of useful functions to obtain information about the physical properties of the system under analysis.
They are the following:
```
nanoribbonS, nanoribbonSA, nanoribbonSZ, Params,  modelS, rectangle_weaklink, rectangle_randombounds_sc, ldosonlattice_averaged_sc
```
Detailed information on any of them can be accessed as follows: `?function_name`


