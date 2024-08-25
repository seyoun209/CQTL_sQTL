# Supp table3 - Pathway and GO terms table
library(data.table)
library(openxlsx)

# Function to write a data frame to a worksheet with formatting
write_formatted_sheet <- function(wb, sheet_name, data) {
  # Add worksheet
  addWorksheet(wb, sheet_name)
  
  # Write data
  writeData(wb, sheet_name, data, startRow = 1, startCol = 1)
  
  # Create styles
  headerStyle <- createStyle(textDecoration = "bold", halign = "center", valign = "center", fgFill = "grey90")
  dataStyle <- createStyle(wrapText = TRUE, halign = "left", valign = "center")

  
  # Apply header style
  addStyle(wb, sheet_name, headerStyle, rows = 1, cols = 1:ncol(data), gridExpand = TRUE)
  
  # Apply data style to all cells except header
  addStyle(wb, sheet_name, dataStyle, rows = 2:(nrow(data)+1), cols = 1:ncol(data), gridExpand = TRUE)

  
  # Set column widths
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
  
  # Freeze the top row
  freezePane(wb, sheet_name, firstRow = TRUE)
}




# 
Go_table <- fread("output/clu_fnf/table/GO_sig.csv")
Pathway_table <- fread("output/clu_fnf/table/pathway_table_sig.csv")



Go_table_OA <- fread("output/clu_oa/table/GO_sig.csv")
Pathway_table_OA <- fread("output/clu_oa/table/pathway_table_sig.csv")

# Create a new workbook
wb <- createWorkbook()
write_formatted_sheet(wb, "Pathway (PBS vs FN-f)", Pathway_table)
write_formatted_sheet(wb, "Pathway (PBS vs OA)", Pathway_table_OA)
write_formatted_sheet(wb, "GO (PBS vs FN-f)", Go_table)
write_formatted_sheet(wb, "GO (PBS vs OA)", Go_table_OA)
saveWorkbook(wb, "output/results_plots/Supp_Data/Supplementary_Data3.xlsx", overwrite = TRUE)
