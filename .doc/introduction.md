### Introduction

This document describes example problems for MODFLOW 6. The examples demonstrate use of the capabilities for select components of MODFLOW 6. Examples have been included for the MODFLOW 6 components summarized in the tables below.

#### Discretization types

| Package | Examples                                                           |
|---------|--------------------------------------------------------------------|
| [DIS](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-dis.html) | ex-gwf-twri01, ex-gwf-lgr, ex-gwf-csub-p01, ex-gwf-csub-p02a, ex-gwf-csub-p02b, ex-gwf-csub-p02c, ex-gwt-keating, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dsupp82, ex-gwt-mt3dms-p01a, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07, ex-gwt-prudic2004t2 |
| [DISV](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-disv.html) | ex-gwt-moc3d-p02tg |


#### Exchanges

| Package | Examples                                                           |
|---------|--------------------------------------------------------------------|
| [GWFGWF](https://modflow6.readthedocs.io/en/latest/_mf6io/exg-gwfgwf.html) | ex-gwf-lgr |
| [GWFGWT](https://modflow6.readthedocs.io/en/latest/_mf6io/exg-gwfgwt.html) | ex-gwt-mt3dms-p01a, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07 |


#### Groundwater Flow Model Internal Flow Packages

| Package | Examples                                                           |
|---------|--------------------------------------------------------------------|
| [NPF](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-npf.html) | ex-gwf-twri01, ex-gwf-lgr, ex-gwf-csub-p01, ex-gwf-csub-p02a, ex-gwf-csub-p02b, ex-gwf-csub-p02c, ex-gwt-keating, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dsupp82, ex-gwt-mt3dms-p01a, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07, ex-gwt-prudic2004t2 |
| [STO](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-sto.html) | ex-gwf-csub-p01, ex-gwf-csub-p02a, ex-gwf-csub-p02b, ex-gwf-csub-p02c, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06 |
| [CSUB](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-csub.html) | ex-gwf-csub-p01, ex-gwf-csub-p02a, ex-gwf-csub-p02b, ex-gwf-csub-p02c |


#### Groundwater Flow Model Standard Boundary Packages

| Package | Examples                                                           |
|---------|--------------------------------------------------------------------|
| [RCH](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-rch.html) | ex-gwt-keating |
| [RCHA](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-rcha.html) | ex-gwf-twri01, ex-gwt-prudic2004t2 |
| [CHD](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-chd.html) | ex-gwf-twri01, ex-gwf-lgr, ex-gwf-csub-p02a, ex-gwf-csub-p02b, ex-gwf-csub-p02c, ex-gwt-keating, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dsupp82, ex-gwt-mt3dms-p01a, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07, ex-gwt-prudic2004t2 |
| [DRN](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-drn.html) | ex-gwf-twri01 |
| [WEL](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-wel.html) | ex-gwf-twri01, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07 |


#### Groundwater Flow Model Advanced Boundary Packages

| Package | Examples                                                           |
|---------|--------------------------------------------------------------------|
| [LAK](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-lak.html) | ex-gwt-prudic2004t2 |
| [SFR](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-sfr.html) | ex-gwf-lgr, ex-gwt-prudic2004t2 |
| [MAW](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-maw.html) | ex-gwt-mt3dsupp82 |
| [MVR](https://modflow6.readthedocs.io/en/latest/_mf6io/gwf-mvr.html) | ex-gwt-mt3dsupp82, ex-gwt-prudic2004t2 |


#### Groundwater Transport Model Internal Flow Packages

| Package | Examples                                                           |
|---------|--------------------------------------------------------------------|
| [ADV](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-adv.html) | ex-gwt-keating, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dsupp82, ex-gwt-mt3dms-p01a, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07, ex-gwt-prudic2004t2 |
| [DSP](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-dsp.html) | ex-gwt-keating, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dsupp82, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07, ex-gwt-prudic2004t2 |
| [MST](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-mst.html) | ex-gwt-keating, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dsupp82, ex-gwt-mt3dms-p01a, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07, ex-gwt-prudic2004t2 |
| [IST](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-ist.html) | ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c |
| [FMI](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-fmi.html) | ex-gwt-keating, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dsupp82, ex-gwt-prudic2004t2 |




#### Groundwater Transport Model Standard Boundary Packages

| Package | Examples                                                           |
|---------|--------------------------------------------------------------------|
| [SSM](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-ssm.html) | ex-gwt-keating, ex-gwt-moc3d-p01a, ex-gwt-moc3d-p01b, ex-gwt-moc3d-p01c, ex-gwt-moc3d-p01d, ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg, ex-gwt-mt3dsupp631, ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dsupp82, ex-gwt-mt3dms-p01a, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p04a, ex-gwt-mt3dms-p04b, ex-gwt-mt3dms-p04c, ex-gwt-mt3dms-p05, ex-gwt-mt3dms-p06, ex-gwt-mt3dms-p07, ex-gwt-prudic2004t2 |
| [SRC](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-src.html) | ex-gwt-moc3d-p02, ex-gwt-moc3d-p02tg |
| [CNC](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-cnc.html) | ex-gwt-mt3dsupp632a, ex-gwt-mt3dsupp632b, ex-gwt-mt3dsupp632c, ex-gwt-mt3dms-p01a, ex-gwt-mt3dms-p01b, ex-gwt-mt3dms-p01c, ex-gwt-mt3dms-p01d, ex-gwt-mt3dms-p03, ex-gwt-mt3dms-p05, ex-gwt-prudic2004t2 |


#### Groundwater Transport Model Advanced Boundary Packages

| Package | Examples                                                           |
|---------|--------------------------------------------------------------------|
| [LKT](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-lkt.html) | ex-gwt-prudic2004t2 |
| [SFT](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-sft.html) | ex-gwt-prudic2004t2 |
| [MWT](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-mwt.html) | ex-gwt-mt3dsupp82 |
| [MVT](https://modflow6.readthedocs.io/en/latest/_mf6io/gwt-mvt.html) | ex-gwt-mt3dsupp82, ex-gwt-prudic2004t2 |


