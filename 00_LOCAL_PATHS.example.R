# LOCAL ONLY: make a private copy named 00_LOCAL_PATHS.R and fill the paths.
# This example file is safe to upload. Never upload the private copy or input data.
Sys.setenv(
  PD_TRPV1_ROOT = "YOUR_PATH_TO_PPMI_RNAseq_IR3_Analysis",
  PD_TRPV1_CLINICAL_DIR = "YOUR_PATH_TO_PPMI_CLINICAL",
  PD_TRPV1_GSE8397_DIR = "YOUR_PATH_TO_GSE8397",
  HOSPITAL_CORRELATION_XLSX = "YOUR_PRIVATE_HOSPITAL_CORRELATION_XLSX"
)
# Optional: Sys.setenv(HOSPITAL_CORRELATION_OUTPUT_DIR="YOUR_PRIVATE_RESULTS_FOLDER")
