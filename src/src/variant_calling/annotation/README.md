# Variant Annotation Module

## Overview

This module provides comprehensive functional and clinical annotation of genomic variants using industry-standard tools and databases. The annotation pipeline integrates multiple annotation sources to deliver rich, clinically-relevant variant information.

### Key Features

- **Functional Annotation**: VEP (Variant Effect Predictor) and SnpEff for consequence prediction
- **Clinical Significance**: ClinVar integration for known pathogenic variants
- **Population Frequencies**: gnomAD database integration for allele frequencies
- **Pathogenicity Prediction**: dbNSFP scores (SIFT, PolyPhen, CADD, etc.)
- **Flexible Output Formats**: VCF, TSV, JSON support
- **Batch Processing**: High-throughput annotation workflows

## Tools and Databases

### Primary Annotation Engines
- **VEP (Variant Effect Predictor)**: Ensembl's comprehensive annotation tool
- **SnpEff**: Fast genetic variant annotation and effect prediction

### Reference Databases
- **dbNSFP**: Database of non-synonymous functional predictions
- **ClinVar**: Clinical significance annotations
- **gnomAD**: Genome/exome population frequency data
- **COSMIC**: Cancer mutation database
- **dbSNP**: Common genetic variants

## Usage Examples

### Basic VEP Annotation

```bash
# Annotate VCF with VEP
python scripts/vep_annotate.py \
  --input variants.vcf \
  --output annotated_variants.vcf \
  --genome GRCh38 \
  --cache-dir /data/vep_cache
```

### SnpEff Annotation

```bash
# Run SnpEff annotation
python scripts/snpeff_annotate.py \
  --input variants.vcf \
  --output snpeff_annotated.vcf \
  --genome hg38 \
  --config config/snpeff.config
```

### Comprehensive Multi-Database Annotation

```bash
# Full annotation pipeline
python scripts/comprehensive_annotate.py \
  --input variants.vcf \
  --output comprehensive_annotated.vcf \
  --vep-cache /data/vep_cache \
  --snpeff-db hg38 \
  --dbnsfp /data/dbNSFP.txt.gz \
  --clinvar /data/clinvar.vcf.gz \
  --gnomad /data/gnomad.vcf.gz \
  --threads 8
```

### Custom Annotation Fields

```bash
# Specify custom annotation fields
python scripts/custom_annotate.py \
  --input variants.vcf \
  --output custom_annotated.vcf \
  --fields "Gene,Consequence,SIFT_score,PolyPhen_score,ClinVar_CLNSIG,gnomAD_AF" \
  --format TSV
```

### Batch Processing

```bash
# Process multiple VCF files
python scripts/batch_annotate.py \
  --input-dir input_vcfs/ \
  --output-dir annotated_vcfs/ \
  --config config/annotation.yaml \
  --parallel 4
```

## Directory Structure

```
annotation/
├── README.md
├── scripts/
│   ├── vep_annotate.py              # VEP annotation wrapper
│   ├── snpeff_annotate.py           # SnpEff annotation wrapper
│   ├── comprehensive_annotate.py    # Multi-tool annotation pipeline
│   ├── custom_annotate.py           # Custom field annotation
│   ├── batch_annotate.py           # Batch processing script
│   ├── filter_annotated.py         # Post-annotation filtering
│   └── utils/
│       ├── annotation_utils.py     # Common annotation functions
│       ├── database_utils.py       # Database interaction utilities
│       └── format_converter.py     # Format conversion utilities
├── config/
│   ├── annotation.yaml             # Main annotation configuration
│   ├── vep.config                  # VEP-specific settings
│   ├── snpeff.config              # SnpEff configuration
│   ├── database_paths.yaml        # Database file paths
│   └── field_mappings.yaml        # Custom field mapping definitions
├── databases/
│   ├── download_databases.sh      # Database download script
│   ├── update_databases.sh        # Database update script
│   └── README_databases.md        # Database setup instructions
├── output/
│   ├── annotated/                 # Annotated VCF files
│   ├── reports/                   # Annotation summary reports
│   └── filtered/                  # Filtered annotation results
├── logs/
│   └── annotation_logs/           # Annotation process logs
└── tests/
    ├── test_annotation.py         # Unit tests
    ├── test_data/                 # Test datasets
    └── integration_tests/         # Integration test suite
```

## Configuration

### Main Configuration File (`config/annotation.yaml`)

```yaml
# Annotation Configuration
genome_build: "GRCh38"
threads: 8
temp_dir: "/tmp/annotation"

# Tool Settings
vep:
  cache_dir: "/data/vep_cache"
  plugins: ["dbNSFP", "CADD", "G2P"]
  fields: ["Gene", "Consequence", "IMPACT", "SYMBOL"]

snpeff:
  database: "hg38"
  config_file: "config/snpeff.config"
  upstream: 5000
  downstream: 5000

# Database Paths
databases:
  dbnsfp: "/data/dbNSFP4.3a.txt.gz"
  clinvar: "/data/clinvar.vcf.gz"
  gnomad: "/data/gnomad.genomes.v3.1.sites.vcf.gz"
  cosmic: "/data/cosmic.vcf.gz"

# Output Settings
output:
  format: "VCF"  # Options: VCF, TSV, JSON
  compress: true
  index: true
```

### Database Setup

```bash
# Download and setup annotation databases
cd databases/
./download_databases.sh

# Update databases (run monthly)
./update_databases.sh
```

## Integration with Other Modules

### Workflow Integration

```python
# Integration with variant calling pipeline
from variant_calling.annotation.scripts.comprehensive_annotate import AnnotationPipeline
from variant_calling.filtering.scripts.variant_filter import VariantFilter

# Annotate called variants
annotator = AnnotationPipeline(config_file="config/annotation.yaml")
annotated_vcf = annotator.annotate("called_variants.vcf")

# Apply post-annotation filters
filter_tool = VariantFilter()
filtered_vcf = filter_tool.filter_by_annotation(annotated_vcf, 
                                               min_qual=20, 
                                               clinvar_pathogenic=True)
```

### Script Chaining

```bash
# Complete variant processing workflow
./variant_calling/scripts/gatk_call_variants.py --input aligned.bam --output raw_variants.vcf
./variant_calling/annotation/scripts/comprehensive_annotate.py --input raw_variants.vcf --output annotated_variants.vcf
./variant_calling/filtering/scripts/variant_filter.py --input annotated_variants.vcf --output final_variants.vcf
./variant_calling/reporting/scripts/generate_report.py --input final_variants.vcf --output variant_report.html
```

## Best Practices

### Performance Optimization

1. **Use Local Caches**: Download and maintain local VEP/SnpEff caches
2. **Parallel Processing**: Utilize multiple threads for large datasets
3. **Selective Annotation**: Only annotate fields needed for downstream analysis
4. **Batch Processing**: Process multiple samples together when possible

### Quality Control

1. **Validate Input**: Check VCF format and chromosome naming conventions
2. **Monitor Resources**: Track memory and disk usage during annotation
3. **Version Tracking**: Document database versions used for reproducibility
4. **Output Validation**: Verify annotation completeness and accuracy

### Database Management

1. **Regular Updates**: Update annotation databases monthly
2. **Version Control**: Maintain multiple database versions for consistency
3. **Backup Strategy**: Keep backups of stable database versions
4. **Access Patterns**: Optimize database access for your typical workload

## Output Formats and Interpretation

### VCF Output

Annotated VCF files include additional INFO fields:

```
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations from SnpEff">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance from ClinVar">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency from gnomAD">
```

### TSV Output

Tab-delimited format with customizable columns:

```
CHROM	POS	REF	ALT	GENE	CONSEQUENCE	SIFT_SCORE	POLYPHEN_SCORE	CLINVAR_CLNSIG	GNOMAD_AF
1	12345	A	T	BRCA1	missense_variant	0.02	0.95	Pathogenic	0.0001
```

## Troubleshooting

### Common Issues

1. **Database Not Found**: Check database paths in configuration
2. **Memory Issues**: Reduce batch size or increase available memory
3. **Slow Performance**: Use local databases and enable parallel processing
4. **Format Errors**: Validate input VCF format and chromosome naming

### Log Analysis

```bash
# Check annotation logs
tail -f logs/annotation_logs/annotation.log

# Search for specific errors
grep -i "error" logs/annotation_logs/*.log
```

## Related Documentation

### Module Documentation
- [Main Pipeline README](../../../README.md) - Overall pipeline documentation
- [Variant Calling README](../README.md) - Variant calling module overview
- [Filtering README](../filtering/README.md) - Post-annotation filtering
- [Quality Control README](../qc/README.md) - Variant quality assessment
- [Reporting README](../reporting/README.md) - Annotation report generation

### External Resources
- [VEP Documentation](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [SnpEff Documentation](http://pcingola.github.io/SnpEff/)
- [ClinVar Documentation](https://www.ncbi.nlm.nih.gov/clinvar/)
- [gnomAD Documentation](https://gnomad.broadinstitute.org/)
- [dbNSFP Documentation](https://sites.google.com/site/jpopgen/dbNSFP)

### Database Resources
- [Database Setup Guide](databases/README_databases.md)
- [VEP Cache Installation](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html)
- [SnpEff Database Download](http://pcingola.github.io/SnpEff/se_buildingdb/)

---

*For technical support or questions about variant annotation, please consult the main pipeline documentation or contact the development team.*
