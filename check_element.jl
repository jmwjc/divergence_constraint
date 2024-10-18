using XLSX

ind = 100
XLSX.openxlsx("./xlsx/check_element.xlsx", mode = "rw") do xf
    sheet = xf[1]
    for i in 1:ind
        for j in 1:ind
    row = "A"
    row_L₂_𝒖 = "B"
    row_Hₑ_𝒖 = "C"
    row_Hₑ_dev = "D"
    row_L₂_𝑝 = "E"
    for (n,L₂_𝒖_,Hₑ_𝒖_,Hₑ_dev_,L₂_𝑝_) in zip(indices,L₂_𝒖,Hₑ_𝒖,Hₑ_dev,L₂_𝑝)
        sheet[row*string(n)] = n
        sheet[row_L₂_𝒖*string(n)] = L₂_𝒖_
        sheet[row_Hₑ_𝒖*string(n)] = Hₑ_𝒖_
        sheet[row_Hₑ_dev*string(n)] = Hₑ_dev_
        sheet[row_L₂_𝑝*string(n)] = L₂_𝑝_
    end
end