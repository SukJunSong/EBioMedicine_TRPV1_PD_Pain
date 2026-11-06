# Input/output paths are configured via environment variables: PD_TRPV1_ROOT.
if (!nzchar(Sys.getenv("PD_TRPV1_ROOT", unset=""))) stop("Set PD_TRPV1_ROOT before running this script (see README.md).", call. = FALSE)

# REVIEWED static audit 2026-11-05; original: 20_Hospital_PBMC_Correlation_Heatmap.R
# Requires original project datasets; full R execution not validated here.
library(readxl)
library(Hmisc)
library(corrplot)

# 데이터 불러오기
# Hospital data must remain local. Set this before running the script, e.g.:
# Set HOSPITAL_CORRELATION_XLSX to a private local Excel file (see README).
input_file <- Sys.getenv("HOSPITAL_CORRELATION_XLSX", unset="")
if (!nzchar(input_file) || !file.exists(input_file)) {
  stop("Set HOSPITAL_CORRELATION_XLSX to an existing local correlation.xlsx; do not publish hospital patient data.")
}
output_dir <- Sys.getenv("HOSPITAL_CORRELATION_OUTPUT_DIR", unset=file.path(Sys.getenv("PD_TRPV1_ROOT", unset=""), "Hospital_PBMC_results"))
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
data <- readxl::read_excel(input_file)

# 숫자형 변수만 선택
data_num <- data[, sapply(data, is.numeric)]

# TRPV1 맨 앞으로
if("TRPV1" %in% colnames(data_num)){
  
  new_order <- c(
    "TRPV1",
    setdiff(colnames(data_num), "TRPV1")
  )
  
  data_num <- data_num[, new_order]
}

# Prism 스타일 Pearson correlation
# pairwise complete observation
cor_res <- rcorr(as.matrix(data_num), type = "pearson")

cor_mat <- cor_res$r
p_mat <- cor_res$P

# significance
sig_level <- matrix("", nrow(p_mat), ncol(p_mat))

sig_level[p_mat < 0.05]  <- "*"
sig_level[p_mat < 0.01]  <- "**"
sig_level[p_mat < 0.001] <- "***"

# 색상
col_fun <- colorRampPalette(c("#648FFF", "white", "#FA8072"))(200)

# PDF 저장
pdf(
  file.path(output_dir, "TRPV1_Correlation_heatmap.pdf"),
  width = 10,
  height = 10
)

# 기존 label 위치 그대로 두고, text 기본 font를 bold로 설정
par(font = 2)

corrplot(
  cor_mat,
  method = "color",
  type = "upper",
  order = "original",
  
  col = col_fun,
  
  addCoef.col = "black",
  number.cex = 0.7,
  number.digits = 2,
  number.font = 1,   # correlation 숫자는 기존처럼 plain 유지
  
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,
  
  diag = FALSE
)

# significance star
for(i in 1:nrow(cor_mat)){
  for(j in 1:ncol(cor_mat)){
    
    if(i < j){
      
      text(
        x = j,
        y = nrow(cor_mat) - i + 1 + 0.22,
        labels = sig_level[i,j],
        cex = 1.1,
        font = 2
      )
      
    }
  }
}

dev.off()
