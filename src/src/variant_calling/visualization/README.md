# Variant Visualization Module

## Overview

This module provides comprehensive visual exploration and analysis capabilities for annotated genomic variants. The visualization pipeline integrates multiple visualization tools and platforms to deliver interactive and publication-ready visualizations of variant data, enabling effective interpretation and communication of genomic findings.

### Key Features

• Interactive Genomic Viewers: IGV and IGV.js integration for detailed variant inspection
• Statistical Dashboards: Web-based dashboards for variant distribution analysis  
• Circos Plots: Genome-wide variant distribution and structural variant visualization
• Manhattan Plots: GWAS-style association and significance visualization
• Heatmaps: Variant density, pathogenicity scores, and annotation correlation matrices
• Multi-sample Comparisons: Comparative visualization across samples and cohorts
• Export Capabilities: High-resolution publication-ready graphics and interactive web formats

## Tools and Platforms

### Primary Visualization Engines

• IGV (Integrative Genomics Viewer): Desktop application for detailed variant inspection
• IGV.js: Web-based genomic visualization for interactive browser-based analysis
• Circos: Circular genome visualization for structural variants and genome-wide patterns
• R/ggplot2: Statistical plotting and publication-quality graphics
• Plotly: Interactive web-based visualizations and dashboards
• D3.js: Custom interactive visualizations

### Dashboard Platforms

• Shiny: R-based interactive web applications
• Dash: Python-based analytical web applications  
• Streamlit: Rapid prototyping of data visualization apps
• Jupyter Notebooks: Exploratory data analysis and visualization

## Usage Examples

### Basic IGV Visualization

```bash
# Launch IGV with variant data
python scripts/igv_visualize.py \
  --input annotated_variants.vcf.gz \
  --reference genome.fa \
  --regions regions_of_interest.bed \
  --output igv_session.xml
```

### IGV.js Web Interface

```bash
# Generate IGV.js web visualization
python scripts/igvjs_generate.py \
  --input annotated_variants.vcf.gz \
  --reference genome.fa \
  --output web_visualization/ \
  --port 8080
```

### Circos Genome Plot

```bash
# Create circos plot for structural variants
python scripts/circos_plot.py \
  --input structural_variants.vcf \
  --config config/circos_config.conf \
  --output circos_plot.png \
  --format png,svg
```

### Manhattan Plot Generation

```bash
# Generate Manhattan plot for variant significance
Rscript scripts/manhattan_plot.R \
  --input gwas_results.tsv \
  --output manhattan_plot.png \
  --title "Genome-wide Variant Associations" \
  --pvalue-column P_VALUE
```

### Variant Density Heatmap

```bash
# Create variant density heatmap
python scripts/variant_heatmap.py \
  --input annotated_variants.vcf.gz \
  --regions chromosomes.bed \
  --output variant_density_heatmap.png \
  --bin-size 1000000
```

### Interactive Dashboard

```bash
# Launch Shiny dashboard
Rscript scripts/launch_dashboard.R \
  --input annotated_variants.vcf.gz \
  --port 3838 \
  --config config/dashboard_config.yaml
```

## Directory Structure

```
visualization/
├── README.md
├── scripts/
│   ├── igv_visualize.py             # IGV session generation
│   ├── igvjs_generate.py            # IGV.js web interface
│   ├── circos_plot.py               # Circos visualization
│   ├── manhattan_plot.R             # Manhattan plot generation
│   ├── variant_heatmap.py           # Heatmap visualization
│   ├── launch_dashboard.R           # Interactive dashboard launcher
│   ├── multi_sample_compare.py      # Multi-sample comparison plots
│   ├── pathogenicity_plots.R        # Pathogenicity score visualizations
│   └── utils/
│       ├── plot_utils.py            # Common plotting functions
│       ├── color_schemes.py         # Standardized color palettes
│       └── export_utils.py          # Export format utilities
├── config/
│   ├── visualization_config.yaml   # Main visualization settings
│   ├── igv_config.xml              # IGV session templates
│   ├── circos_config.conf          # Circos plot configuration
│   ├── dashboard_config.yaml       # Dashboard layout settings
│   └── color_palettes.yaml        # Color scheme definitions
├── templates/
│   ├── igv_session_template.xml    # IGV session template
│   ├── dashboard_template.html     # Dashboard HTML template
│   ├── report_template.Rmd        # R Markdown report template
│   └── circos_template.conf       # Circos configuration template
├── web/
│   ├── dashboard/                  # Web dashboard files
│   │   ├── app.R                  # Shiny application
│   │   ├── ui.R                   # User interface definition
│   │   ├── server.R               # Server logic
│   │   └── www/                   # Static web assets
│   └── igvjs/                     # IGV.js web viewer
│       ├── index.html             # Main IGV.js page
│       ├── js/                    # JavaScript files
│       └── css/                   # Styling files
├── output/
│   ├── plots/                     # Generated static plots
│   ├── interactive/               # Interactive visualizations
│   ├── sessions/                  # IGV session files
│   └── reports/                   # Visualization reports
├── data/
│   ├── reference/                 # Reference data for visualization
│   ├── annotations/               # Additional annotation tracks
│   └── examples/                  # Example datasets
└── tests/
    ├── test_visualization.py      # Unit tests
    ├── test_data/                 # Test datasets
    └── integration_tests/         # Integration test suite
```

## Configuration

### Main Configuration File (config/visualization_config.yaml)

```yaml
# Visualization Configuration

# General Settings
genome_build: "GRCh38"
output_format: ["png", "svg", "pdf", "html"]
resolution: 300  # DPI for static plots
interactive: true

# IGV Settings
igv:
  memory: "4g"
  port: 60151
  batch_mode: true
  snapshot_format: "png"
  
# Circos Settings  
circos:
  image_size: "3000x3000"
  background_color: "white"
  highlight_variants: true
  
# Dashboard Settings
dashboard:
  port: 3838
  theme: "bootstrap"
  max_file_size: "100MB"
  cache_enabled: true
  
# Plot Settings
plots:
  width: 12
  height: 8  
  dpi: 300
  font_family: "Arial"
  color_palette: "viridis"
  
# Export Settings
export:
  formats: ["png", "pdf", "svg", "html"]
  high_res: true
  vector_graphics: true
```

### IGV Configuration Template

```xml
<?xml version="1.0" encoding="UTF-8"?>
<Session genome="hg38" version="8">
  <Resources>
    <Resource path="annotated_variants.vcf.gz"/>
    <Resource path="reference.fa"/>
  </Resources>
  <Panel name="Panel1">
    <Track name="Variants" visible="true" height="60"/>
    <Track name="Genes" visible="true" height="35"/>
  </Panel>
</Session>
```

## Integration with Annotation Module

### Data Flow Integration

```python
# Integration with annotation pipeline
from variant_calling.annotation.scripts.comprehensive_annotate import AnnotationPipeline
from variant_calling.visualization.scripts.igv_visualize import IGVVisualizer

# Annotate variants and visualize
annotator = AnnotationPipeline(config_file="config/annotation.yaml")
annotated_vcf = annotator.annotate("variants.vcf")

# Generate visualizations
visualizer = IGVVisualizer()
visualization_files = visualizer.create_session(annotated_vcf, 
                                               regions_of_interest="roi.bed",
                                               output_dir="visualizations/")
```

### Workflow Integration

```bash
# Complete annotation to visualization workflow
./annotation/scripts/comprehensive_annotate.py --input variants.vcf --output annotated.vcf
./visualization/scripts/igv_visualize.py --input annotated.vcf --output igv_session.xml
./visualization/scripts/launch_dashboard.R --input annotated.vcf --port 3838
```

## Best Practices

### Performance Optimization

1. Data Preprocessing: Index VCF files and use tabix for efficient random access
2. Efficient Rendering: Use appropriate plot dimensions and resolution for intended use
3. Caching: Cache intermediate visualization data for repeated analyses
4. Memory Management: Process large datasets in chunks for memory-efficient visualization

### Visual Design Guidelines

1. Color Consistency: Use standardized color palettes across all visualizations
2. Accessibility: Ensure color schemes work for colorblind users
3. Information Density: Balance detail with readability in complex plots
4. Interactive Elements: Provide tooltips and zoom capabilities for detailed exploration

### Data Integration

1. Standardized Formats: Ensure consistent data format expectations across visualization tools
2. Annotation Completeness: Verify all required annotation fields are present before visualization
3. Quality Filtering: Apply appropriate quality filters before visualization to reduce noise
4. Multi-scale Views: Provide both genome-wide and detailed local views

## Visualization Types and Use Cases

### IGV/IGV.js Browser Views

Detailed inspection of individual variants with:
- Read alignment support
- Annotation track overlay  
- Splice site visualization
- Structural variant breakpoints

### Statistical Distribution Plots

Population-level variant analysis:
- Allele frequency distributions
- Quality score distributions
- Variant type proportions
- Sample-level variant counts

### Genome-wide Visualizations

Chromosomal and genome-wide patterns:
- Circos plots for structural variants
- Manhattan plots for association studies
- Karyotype plots with variant density
- Comparative genomics visualizations

### Correlation and Comparison Matrices

Multi-dimensional analysis:
- Sample-to-sample variant overlap heatmaps
- Annotation correlation matrices
- Pathogenicity score comparisons
- Caller concordance visualizations

## Output Formats and Export Options

### Static Graphics
- PNG: Web display and presentations
- PDF: Publication-quality vector graphics
- SVG: Scalable vector graphics for web
- TIFF: High-resolution scientific publications

### Interactive Formats
- HTML: Standalone interactive plots
- Shiny Apps: Full-featured web applications
- Jupyter Notebooks: Exploratory analysis documents
- IGV Sessions: Shareable genomic browser sessions

## Troubleshooting

### Common Issues

1. Memory Issues: Reduce dataset size or increase available memory
2. Slow Rendering: Optimize plot complexity or use data sampling
3. Format Compatibility: Verify input file formats and indices
4. Port Conflicts: Check for conflicting services on dashboard ports

### Performance Tuning

```bash
# Monitor resource usage
./scripts/monitor_visualization.sh

# Optimize large datasets
python scripts/optimize_vcf.py --input large_variants.vcf --output optimized.vcf
```

## Related Documentation

### Module Documentation

• [Main Pipeline README](../../README.md) - Overall pipeline documentation
• [Variant Calling README](../README.md) - Variant calling module overview  
• [Annotation README](../annotation/README.md) - Variant annotation module
• [Filtering README](../filtering/README.md) - Post-annotation filtering
• [Quality Control README](../qc/README.md) - Variant quality assessment
• [Reporting README](../reporting/README.md) - Automated report generation

### External Resources

• [IGV Documentation](https://software.broadinstitute.org/software/igv/)
• [IGV.js Documentation](https://github.com/igvteam/igv.js/)
• [Circos Documentation](http://circos.ca/documentation/)
• [Shiny Documentation](https://shiny.rstudio.com/)
• [Plotly Documentation](https://plotly.com/r/)
• [ggplot2 Documentation](https://ggplot2.tidyverse.org/)

### Visualization Resources

• [Color Brewer](https://colorbrewer2.org/) - Color schemes for maps and charts
• [Data Visualization Guidelines](https://datavisualization.ch/) - Best practices
• [Genomic Visualization Tools](https://academic.oup.com/bioinformatics) - Literature reviews

For technical support or questions about variant visualization, please consult the main pipeline documentation or contact the development team.
