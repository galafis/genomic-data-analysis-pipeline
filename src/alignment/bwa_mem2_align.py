#!/usr/bin/env python3
"""
BWA-MEM2 Alignment Wrapper

Author: Gabriel Demetrios Lafis
Description: Professional wrapper for BWA-MEM2 alignment with comprehensive
             configuration, logging, and quality control integration.
Usage:
    python bwa_mem2_align.py --config config.yaml
    python bwa_mem2_align.py --input R1.fq.gz R2.fq.gz --reference genome.fa --output sample.bam
"""

import argparse
import subprocess
import logging
import sys
import os
from pathlib import Path
from typing import Optional, Tuple
import yaml
import json
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('bwa_mem2_alignment.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)


class BWAMem2Aligner:
    """Professional BWA-MEM2 alignment wrapper with comprehensive features."""
    
    def __init__(self, config: Optional[dict] = None):
        """
        Initialize BWA-MEM2 aligner.
        
        Args:
            config: Configuration dictionary with alignment parameters
        """
        self.config = config or self._default_config()
        self.validate_dependencies()
        
    def _default_config(self) -> dict:
        """Return default configuration parameters."""
        return {
            'threads': 16,
            'mark_duplicates': True,
            'min_mapq': 20,
            'sort_memory': '4G',
            'platform': 'ILLUMINA',
            'create_index': True,
            'collect_metrics': True
        }
    
    def validate_dependencies(self) -> None:
        """Validate that required tools are available."""
        required_tools = ['bwa-mem2', 'samtools']
        for tool in required_tools:
            if subprocess.call(['which', tool], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
                raise RuntimeError(f"Required tool '{tool}' not found in PATH")
        logger.info("All required dependencies validated")
    
    def build_read_group(self, sample_id: str, library: str = 'lib1', platform: str = 'ILLUMINA') -> str:
        """Build read group string for BAM file."""
        rg = f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tLB:{library}\\tPL:{platform}"
        return rg
    
    def run_alignment(self,
                     read1: str,
                     read2: Optional[str],
                     reference: str,
                     output_bam: str,
                     sample_id: str,
                     library: str = 'lib1') -> Tuple[bool, str]:
        """
        Run BWA-MEM2 alignment.
        
        Args:
            read1: Path to R1 FASTQ file
            read2: Path to R2 FASTQ file (None for single-end)
            reference: Path to reference genome
            output_bam: Path to output BAM file
            sample_id: Sample identifier
            library: Library name
            
        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            logger.info(f"Starting alignment for sample: {sample_id}")
            logger.info(f"Reference: {reference}")
            logger.info(f"R1: {read1}")
            if read2:
                logger.info(f"R2: {read2}")
            
            # Build read group
            rg = self.build_read_group(sample_id, library, self.config['platform'])
            
            # Build BWA-MEM2 command
            bwa_cmd = [
                'bwa-mem2', 'mem',
                '-t', str(self.config['threads']),
                '-R', rg,
                reference,
                read1
            ]
            
            if read2:
                bwa_cmd.append(read2)
            
            # Build samtools sort command
            sort_cmd = [
                'samtools', 'sort',
                '-@', str(max(1, self.config['threads'] // 2)),
                '-m', self.config['sort_memory'],
                '-o', output_bam,
                '-'
            ]
            
            logger.info(f"BWA-MEM2 command: {' '.join(bwa_cmd)}")
            logger.info(f"samtools sort command: {' '.join(sort_cmd)}")
            
            # Run alignment pipeline
            with subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as bwa_proc:
                with subprocess.Popen(sort_cmd, stdin=bwa_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as sort_proc:
                    bwa_proc.stdout.close()
                    stdout, stderr = sort_proc.communicate()
                    
                    if sort_proc.returncode != 0:
                        error_msg = f"Alignment failed: {stderr.decode()}" 
                        logger.error(error_msg)
                        return False, error_msg
            
            # Create index if requested
            if self.config['create_index']:
                self.index_bam(output_bam)
            
            # Collect metrics if requested
            if self.config['collect_metrics']:
                self.collect_alignment_metrics(output_bam, reference)
            
            logger.info(f"Alignment completed successfully: {output_bam}")
            return True, f"Alignment successful: {output_bam}"
            
        except Exception as e:
            error_msg = f"Error during alignment: {str(e)}"
            logger.error(error_msg)
            return False, error_msg
    
    def index_bam(self, bam_file: str) -> None:
        """Create BAM index file."""
        logger.info(f"Indexing BAM file: {bam_file}")
        cmd = ['samtools', 'index', '-@', str(self.config['threads']), bam_file]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Indexing failed: {result.stderr}")
            raise RuntimeError(f"Failed to index BAM: {result.stderr}")
        logger.info(f"BAM indexed successfully: {bam_file}.bai")
    
    def collect_alignment_metrics(self, bam_file: str, reference: str) -> None:
        """Collect alignment statistics using samtools."""
        logger.info(f"Collecting alignment metrics for: {bam_file}")
        
        # Get output directory
        out_dir = Path(bam_file).parent / 'qc'
        out_dir.mkdir(parents=True, exist_ok=True)
        
        base_name = Path(bam_file).stem
        
        # Collect flagstat
        flagstat_file = out_dir / f"{base_name}.flagstat"
        cmd = ['samtools', 'flagstat', bam_file]
        with open(flagstat_file, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        logger.info(f"Flagstat saved to: {flagstat_file}")
        
        # Collect stats
        stats_file = out_dir / f"{base_name}.stats"
        cmd = ['samtools', 'stats', bam_file]
        with open(stats_file, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True)
        logger.info(f"Stats saved to: {stats_file}")
        
    def batch_align(self, samples_file: str, reference: str, output_dir: str) -> dict:
        """
        Batch align multiple samples from a samples file.
        
        Args:
            samples_file: TSV file with columns: sample_id, read1, read2, library
            reference: Path to reference genome
            output_dir: Output directory for BAM files
            
        Returns:
            Dictionary with alignment results
        """
        results = {'success': [], 'failed': []}
        
        # Create output directory
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Read samples file
        with open(samples_file, 'r') as f:
            header = f.readline().strip().split('\t')
            for line in f:
                fields = line.strip().split('\t')
                sample_data = dict(zip(header, fields))
                
                sample_id = sample_data['sample_id']
                read1 = sample_data['read1']
                read2 = sample_data.get('read2', None)
                library = sample_data.get('library', 'lib1')
                
                output_bam = os.path.join(output_dir, f"{sample_id}.sorted.bam")
                
                success, msg = self.run_alignment(read1, read2, reference, output_bam, sample_id, library)
                
                if success:
                    results['success'].append(sample_id)
                else:
                    results['failed'].append({'sample': sample_id, 'error': msg})
        
        # Log summary
        logger.info(f"Batch alignment complete: {len(results['success'])} succeeded, {len(results['failed'])} failed")
        return results


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Professional BWA-MEM2 alignment wrapper',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Single sample alignment (paired-end)
  python bwa_mem2_align.py --r1 sample_R1.fq.gz --r2 sample_R2.fq.gz \
    --reference genome.fa --output sample.bam --sample-id sample1
  
  # Batch alignment from samples file
  python bwa_mem2_align.py --batch samples.tsv --reference genome.fa \
    --output-dir aligned/
  
  # With custom configuration
  python bwa_mem2_align.py --config alignment_config.yaml \
    --batch samples.tsv --reference genome.fa --output-dir aligned/
        '''
    )
    
    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--r1', help='R1 FASTQ file (paired-end or single-end)')
    input_group.add_argument('--batch', help='TSV file with sample information')
    
    parser.add_argument('--r2', help='R2 FASTQ file (for paired-end)')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA')
    parser.add_argument('--sample-id', help='Sample identifier (required for single sample)')
    parser.add_argument('--library', default='lib1', help='Library name (default: lib1)')
    
    # Output options
    parser.add_argument('--output', help='Output BAM file (for single sample)')
    parser.add_argument('--output-dir', help='Output directory (for batch mode)')
    
    # Configuration
    parser.add_argument('--config', help='YAML configuration file')
    parser.add_argument('--threads', type=int, help='Number of threads')
    parser.add_argument('--platform', default='ILLUMINA', help='Sequencing platform')
    
    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_arguments()
    
    # Load configuration
    config = None
    if args.config:
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    else:
        config = {}
    
    # Override config with command line arguments
    if args.threads:
        config['threads'] = args.threads
    if args.platform:
        config['platform'] = args.platform
    
    # Initialize aligner
    aligner = BWAMem2Aligner(config)
    
    try:
        # Batch mode
        if args.batch:
            if not args.output_dir:
                raise ValueError("--output-dir required for batch mode")
            
            logger.info(f"Starting batch alignment from: {args.batch}")
            results = aligner.batch_align(args.batch, args.reference, args.output_dir)
            
            # Save results
            results_file = os.path.join(args.output_dir, 'alignment_results.json')
            with open(results_file, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to: {results_file}")
            
            # Exit with error code if any samples failed
            if results['failed']:
                sys.exit(1)
        
        # Single sample mode
        else:
            if not args.sample_id:
                raise ValueError("--sample-id required for single sample mode")
            if not args.output:
                raise ValueError("--output required for single sample mode")
            
            success, msg = aligner.run_alignment(
                args.r1,
                args.r2,
                args.reference,
                args.output,
                args.sample_id,
                args.library
            )
            
            if not success:
                logger.error(msg)
                sys.exit(1)
            
            logger.info("Alignment pipeline completed successfully")
    
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        sys.exit(1)


if __name__ == '__main__':
    main()
