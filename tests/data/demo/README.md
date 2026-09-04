# BioMA public demo data

This directory is deliberately independent of the production cluster.  It is
small enough for CI and for a first installation check, and contains no
`/usr_storage` paths or species-private files.

Create a disposable demo project with:

```bash
python tests/data/demo/make_demo.py --outdir /tmp/bioma-demo
bin/bioma project /tmp/bioma-demo/project.ini --dry-run
```

The default generator writes lightweight raster placeholders.  They are enough
for configuration validation and the module smoke tests.  To make readable
GeoTIFFs for a local raster/GDAL test, install the locked environment and run:

```bash
python tests/data/demo/make_demo.py --outdir /tmp/bioma-demo --valid-rasters
```

The generated project keeps the three genomic streams explicit:

- `adaptive_sites.vcf`: adaptive loci for GF/RONA;
- `whole_genome.mis0.9.maf0.00001.vcf`: filtered whole-genome loci for MAR;
- `load_3vcf/`: four derived/SIFT load VCFs for loadM/loadD.

The 19 `LD_BIO*.prune.in` files are small deterministic fixtures.  They model
the lists normally produced by PLINK and are recorded by the RONA manifest;
they are not inferred from the population environment table.
