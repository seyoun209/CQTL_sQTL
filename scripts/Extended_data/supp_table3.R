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




