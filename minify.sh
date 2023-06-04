#!/bin/sh

cargo equip --lib --mod-path crate::__fps --remove docs --remove comments --minify libs > minified.rs