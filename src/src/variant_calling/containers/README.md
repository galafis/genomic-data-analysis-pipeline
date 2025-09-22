# 🧪 Containers — Variant Calling

Este diretório contém definições de containers para o módulo Variant Calling, garantindo reprodutibilidade, portabilidade e isolamento de dependências em ambientes locais, HPC e nuvem. O padrão visual e organizacional segue os demais diretórios do módulo.

## 📌 Finalidade
- Reprodutibilidade: fixar versões de ferramentas (GATK, FreeBayes, bcftools, VEP, SnpEff) e bases (VEP cache/SnpEff DB).
- Portabilidade: executar scripts e workflows de forma consistente em Docker/Podman e Singularity/Apptainer.
- Isolamento: reduzir conflitos de bibliotecas do sistema e simplificar setup.

## 🗂️ Estrutura
- Dockerfile.gatk — imagem com GATK e ferramentas de suporte.
- Dockerfile.freebayes — imagem para FreeBayes + bcftools/tabix.
- Dockerfile.annotation — imagem com VEP ou SnpEff e utilitários.
- README.md — este arquivo.

## 🧰 Exemplos de Dockerfile

### Dockerfile.gatk (exemplo)
```Dockerfile
FROM ubuntu:22.04
ARG GATK_VERSION=4.4.0.0
RUN apt-get update && apt-get install -y \
    openjdk-11-jre-headless wget curl git python3 python3-pip \
    bcftools tabix && rm -rf /var/lib/apt/lists/*
RUN wget -q https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip \
 && apt-get update && apt-get install -y unzip \
 && unzip gatk-${GATK_VERSION}.zip -d /opt \
 && ln -s /opt/gatk-${GATK_VERSION}/gatk /usr/local/bin/gatk \
 && rm -f gatk-${GATK_VERSION}.zip
ENV JAVA_OPTS="-Xmx8g"
WORKDIR /work
```

### Dockerfile.annotation (VEP) (exemplo)
```Dockerfile
FROM ubuntu:22.04
RUN apt-get update && apt-get install -y \
    perl build-essential cpanminus curl wget git \
    tabix bcftools zlib1g-dev liblzma-dev libbz2-dev \
    python3 python3-pip && rm -rf /var/lib/apt/lists/*
# Instalar VEP (versão fixa)
ARG VEP_VERSION=110
RUN cpanm Archive::Zip Module::Build DBI DBD::mysql JSON \
 && git clone --branch release/${VEP_VERSION} https://github.com/Ensembl/ensembl-vep.git /opt/vep \
 && perl /opt/vep/INSTALL.pl --NO_HTSLIB --AUTO a --PLUGINS all --QUIET
ENV PATH="/opt/vep:/opt/vep/bin:$PATH"
WORKDIR /work
```

## 📦 Build das Imagens
- Docker/Podman:
  - GATK: docker build -t pipeline/gatk:4.4 -f Dockerfile.gatk .
  - FreeBayes: docker build -t pipeline/freebayes:1.3 -f Dockerfile.freebayes .
  - Annotation: docker build -t pipeline/annotation:vep110 -f Dockerfile.annotation .

- Singularity/Apptainer (a partir do Docker):
  - apptainer build gatk_4.4.sif docker-daemon://pipeline/gatk:4.4
  - apptainer build annotation_vep110.sif docker-daemon://pipeline/annotation:vep110

- Singularity a partir de recipe (exemplo Singularity.def):
```def
Bootstrap: docker
From: ubuntu:22.04
%post
  apt-get update && apt-get install -y bcftools tabix && rm -rf /var/lib/apt/lists/*
%runscript
  exec "$@"
```
Build: apptainer build tools.sif Singularity.def

## ▶️ Uso dos Containers
- Docker (montando diretórios de dados):
  - docker run --rm -v $PWD:/work -w /work pipeline/gatk:4.4 \
    gatk HaplotypeCaller -I input.bam -R reference.fa -O output.vcf.gz

- Singularity/Apptainer:
  - apptainer exec -B $PWD:/work tools.sif \
    gatk HaplotypeCaller -I input.bam -R reference.fa -O output.vcf.gz

- Variáveis comuns:
  - VEP_CACHE=/path/to/vep_cache (montar com -v ou -B)
  - SNPEFF_DB=/path/to/snpeff/data

## 🔗 Integração com Scripts e Workflows
- Scripts (scripts/): defina CONTAINER_ENGINE=(docker|podman|apptainer) e IMAGES no topo do script.
  - Ex.: ./scripts/run_gatk_haplotypecaller.sh usa docker run ... se CONTAINER_ENGINE=docker.
- Workflows (workflows/):
  - Nextflow: use executor.container = 'pipeline/gatk:4.4' ou profiles para cada processo.
  - Snakemake: rule X: container: "docker://pipeline/gatk:4.4".
  - CWL: DockerRequirement com imagePull.

## 🔒 Boas Práticas de Segurança e Versionamento
- Execute como usuário não-root quando possível (USER 1000:1000; ajuste permissões).
- Fixe versões/sha256 das imagens base e ferramentas (ARG + checksums; FROM ubuntu@sha256:...).
- Não inclua dados sensíveis nas imagens; monte via volumes.
- Varredura de vulnerabilidades: trivy/grype nas imagens antes do release.
- Mantenha imagens mínimas (limpe caches, use --no-install-recommends, multi-stage quando aplicável).
- Assinatura e proveniência: cosign/slsa para publicar digests assinados.
- Versionamento semântico das tags: ex. pipeline/gatk:4.4.0-1.
- Documente hashes/manifest em containers/MANIFEST.md com imagem → digest.

## ✅ Dicas de Uso
- Monte caches (VEP/SnpEff) como somente leitura quando possível.
- Para HPC sem Docker, prefira Apptainer com --contain e --cleanenv.
- Ajuste threads com variáveis (ex.: GATK --native-pair-hmm-threads)
- Garanta tabix indices após gerar VCFs.

## 🔗 Links Relacionados
- README do módulo: src/src/variant_calling/README.md
- Integração downstream: src/src/variant_calling/integration/README.md
- Workflows: src/src/variant_calling/workflows/README.md
- Scripts: src/src/variant_calling/scripts/README.md
- Pipeline principal: README.md na raiz do repositório

## 📄 Licença
Segue a licença MIT do projeto (LICENSE na raiz).
