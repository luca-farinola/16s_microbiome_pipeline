# Makefile for 16S Pathoscope Pipeline

.PHONY: help setup databases submit monitor clean preview logs

# Default target
help:
	@echo "16S Pathoscope Microbiome Pipeline - Makefile"
	@echo "=============================================="
	@echo ""
	@echo "Available targets:"
	@echo "  make setup              - Initial setup (check tools, create dirs)"
	@echo "  make databases          - Download and prepare SILVA/Pathoscope databases"
	@echo "  make preview            - Preview pipeline (dry-run)"
	@echo "  make submit             - Submit pipeline to SLURM"
	@echo "  make monitor            - Monitor running job"
	@echo "  make logs               - Show pipeline logs"
	@echo "  make clean              - Remove work directory and cache"
	@echo "  make clean-all          - Remove all analysis outputs (careful!)"
	@echo "  make report             - Open execution report"
	@echo ""
	@echo "Configuration:"
	@echo "  Edit nextflow.config to customize database paths"
	@echo "  Edit submit_pipeline.sh to customize SLURM settings"
	@echo ""

# Setup initial environment
setup:
	@echo "Checking dependencies..."
	@command -v nextflow >/dev/null 2>&1 || { echo "Nextflow not found. Installing..."; curl -s https://get.nextflow.io | bash; }
	@command -v slurmctld >/dev/null 2>&1 || echo "Warning: SLURM not found (required for submission)"
	@command -v singularity >/dev/null 2>&1 || echo "Warning: Singularity not found (required for containers)"
	@chmod +x submit_pipeline.sh setup_databases.sh
	@echo "Setup complete!"

# Download and prepare databases
databases:
	@echo "Starting database setup..."
	@chmod +x setup_databases.sh
	@./setup_databases.sh .
	@echo ""
	@echo "Database setup complete! Update nextflow.config with the paths provided."

# Dry run - preview pipeline
preview:
	@echo "Running pipeline preview (dry-run)..."
	@nextflow run main_pathscope_16s.nf -preview

# Submit to SLURM
submit:
	@echo "Submitting pipeline to SLURM..."
	@chmod +x submit_pipeline.sh
	@sbatch submit_pipeline.sh
	@echo ""
	@echo "Job submitted! Check status with: make monitor"

# Monitor running job
monitor:
	@echo "Checking job status..."
	@squeue -u $$USER | grep -E "microbiome|pathoscope" || echo "No running jobs found"
	@echo ""
	@echo "Recent nextflow log:"
	@tail -20 .nextflow.log 2>/dev/null || echo ".nextflow.log not found"

# Show logs
logs:
	@echo "=== Nextflow Log ==="
	@tail -50 .nextflow.log 2>/dev/null || echo "No nextflow log found"
	@echo ""
	@echo "=== SLURM Errors ==="
	@ls -t slurm_*.err 2>/dev/null | head -1 | xargs cat || echo "No SLURM errors found"

# Clean work directory
clean:
	@echo "Cleaning nextflow work directory..."
	@nextflow clean -f
	@echo "Done!"

# Clean everything (careful!)
clean-all: clean
	@echo "WARNING: This will remove all analysis outputs!"
	@read -p "Are you sure? (y/N) " -n 1 -r; \
	echo; \
	if [[ $$REPLY =~ ^[Yy]$$ ]]; then \
		rm -rf ../02.Analysis && echo "Analysis directory removed"; \
	else \
		echo "Cancelled"; \
	fi

# Open execution report
report:
	@echo "Looking for execution reports..."
	@ls -t reports/execution_*.html 2>/dev/null | head -1 | xargs -I {} echo "Opening: {}" && open {}

# Install Nextflow locally
install-nextflow:
	@echo "Installing Nextflow..."
	@curl -s https://get.nextflow.io | bash
	@chmod +x nextflow
	@echo "Nextflow installed successfully!"

# Check configuration
check-config:
	@echo "Checking configuration..."
	@echo ""
	@echo "Nextflow version:"
	@nextflow -version
	@echo ""
	@echo "Raw data directory:"
	@ls -ld ../01.RawData 2>/dev/null || echo "Not found!"
	@echo ""
	@echo "Sample count:"
	@find ../01.RawData -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l
	@echo ""
	@echo "Database paths in config:"
	@grep -E "silva_db|pathscope_ref" nextflow.config || echo "Not configured!"

# Validate inputs
validate-inputs:
	@echo "Validating input data..."
	@echo ""
	@echo "Checking for paired FASTQ files..."
	@find ../01.RawData -name "*R1*.fastq.gz" 2>/dev/null | wc -l | xargs echo "Found R1 files:"
	@find ../01.RawData -name "*R2*.fastq.gz" 2>/dev/null | wc -l | xargs echo "Found R2 files:"
	@echo ""
	@echo "Checking file integrity..."
	@for f in ../01.RawData/*/*fastq.gz; do \
		if ! gzip -t $$f 2>/dev/null; then echo "CORRUPT: $$f"; fi; \
	done && echo "All files OK"

# Quick start
quick-start:
	@echo "Quick Start for 16S Pathoscope Pipeline"
	@echo "========================================"
	@echo ""
	@echo "1. Setting up environment..."
	@make setup
	@echo ""
	@echo "2. Preparing databases..."
	@echo "   Run: make databases"
	@echo ""
	@echo "3. Configure paths..."
	@echo "   Edit: nextflow.config"
	@echo "   Edit: submit_pipeline.sh"
	@echo ""
	@echo "4. Submit job..."
	@echo "   Run: make submit"
	@echo ""
	@echo "5. Monitor progress..."
	@echo "   Run: make monitor"
	@echo ""
	@echo "Documentation: README.md, QUICKSTART.md"

# View pipeline structure
info:
	@echo "16S Pathoscope Pipeline Structure"
	@echo "=================================="
	@echo ""
	@echo "Work Directory: work/"
	@ls -la work 2>/dev/null | tail -5 || echo "Not created yet"
	@echo ""
	@echo "Configuration files:"
	@ls -1 *.config *.yml 2>/dev/null || echo "No config files found"
	@echo ""
	@echo "Pipeline scripts:"
	@ls -1 *.nf 2>/dev/null || echo "No Nextflow scripts found"
	@echo ""
	@echo "Documentation:"
	@ls -1 *.md 2>/dev/null | head -5

.DEFAULT_GOAL := help

