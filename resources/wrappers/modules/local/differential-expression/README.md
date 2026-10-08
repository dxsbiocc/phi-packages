# Differential-expression wrapper dependencies

The eight wrappers retain their analysis scripts, parameter contracts and entrypoints.
Every module declares its executable R dependencies in its own `environment.yml`;
Nextflow's `conda` profile resolves that module-local file.

Docker and Singularity/Apptainer use a frozen public dependency image containing the
union of those module environments. No locally built `phi/` image or image recipe
directory is required. The public dependency image is provisioned from Conda packages,
without uploading analysis scripts or project data.

## Fixed container references

- Docker/OCI: `community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned@sha256:e500cffd28caf431b76344610bd213a30b536e2d2fa58418255ec00238e3ab12`
- Singularity/Apptainer: `oras://community.wave.seqera.io/library/bioconductor-biocparallel_bioconductor-deseq2_bioconductor-edger_bioconductor-limma_pruned:ccbe9c69fe5b6749`
- Platform: Linux AMD64. On ARM, the OCI image requires Docker's AMD64 emulation;
  the Conda profile resolves native packages when available. Native SIF pulls require
  a Singularity/Apptainer version with ORAS support.
- Image build input environment SHA-256: `94bcfa837334118dd287c271aa163087f416a56547c59cf42f691eeb79930cb3`.

R and Bioconductor version pins retain the compatible cohort recorded by the previous
runtime. DESeq2 explicitly includes `r-ashr` for its existing shrinkage method, and
DREAM explicitly includes the directly used `r-lme4` dependency. Both versions match
the original local runtime (`ashr` 2.2.63 and `lme4` 2.0.6) and are pinned in the
module-local specifications.

The Docker image build loads all nine required R namespaces in its final runtime stage;
its successful build validates the installed libraries. The fixed OCI reference uses the
verified manifest digest. The image uses package-sized layers to keep downloads bounded.

The RNA-seq workflow continues to include this family's QC module. Its runtime profiles
are unchanged; the QC module itself provides the new container and Conda references.

## Updating dependencies

After changing module dependencies, provision both frozen formats from their union,
wait for successful builds, check required R libraries, and replace the fixed references
in all eight process declarations together. Do not copy Phi's core runtime environment
catalog into this content repository.

Provisioning documentation: [Wave packages](https://docs.seqera.io/wave/provisioning),
[API reference](https://docs.seqera.io/wave/api),
[public community images](https://docs.seqera.io/wave/seqera-containers).
