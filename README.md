# SV-NGS-Nextflow

**SV-NGS-Nextflow** is a Nextflow-based pipeline designed for the detection and analysis of structural variants (SVs) from next-generation sequencing (NGS) data.

## Features

- **Alignment**: Supports BWA and Minimap2 for read alignment.
- **Variant Calling**: Utilizes Manta and SvisionPro for SV detection.
- **Annotation**: Employs SnpEff for variant annotation.
- **Comparative Analysis**: Compares results between Manta and SvisionPro.

## Requirements

- **Nextflow**: Ensure Nextflow is installed. [Installation Guide](https://www.nextflow.io/docs/latest/getstarted.html)
- **Java**: Nextflow requires Java 8 or higher.
- **Docker**: (Optional) For containerized execution.

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/DarwinianCoder/SV-NGS-Nextflow-.git
   cd SV-NGS-Nextflow-

