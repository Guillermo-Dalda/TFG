$output_file = "output.txt"

Clear-Content -Path $output_file

$final_output = ""

for ($i = 1; $i -le 30; $i++) {
    $row = & $args[0]
    $row = $row -replace '\t+', "`t" -replace '\n', '' -replace '\.', ','

    $final_output += $row + "`n"
}
Add-Content -Path $output_file -Value $final_output