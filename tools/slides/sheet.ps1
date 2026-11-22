# A labelled contact sheet of every PNG in a folder (#1071: picking slideshow shots). Usage:
#   powershell -File tools/slides/sheet.ps1 -Dir <folder> -Out <sheet.png> [-Across 4] [-Width 480]
param([string]$Dir, [string]$Out, [int]$Across = 4, [int]$Width = 480)
Add-Type -AssemblyName System.Drawing
$files = Get-ChildItem -Path $Dir -Filter *.png | Where-Object { $_.FullName -ne (Resolve-Path -LiteralPath $Dir).Path + '\' + (Split-Path $Out -Leaf) } | Sort-Object Name
if ($files.Count -eq 0) { Write-Output "sheet: no PNGs in $Dir"; exit 1 }
$first = [System.Drawing.Image]::FromFile($files[0].FullName)
$h = [int]($Width * $first.Height / $first.Width); $first.Dispose()
$label = 22
$rows = [math]::Ceiling($files.Count / $Across)
$sheet = New-Object System.Drawing.Bitmap ($Across * $Width), ($rows * ($h + $label))
$g = [System.Drawing.Graphics]::FromImage($sheet)
$g.Clear([System.Drawing.Color]::FromArgb(24, 24, 24))
$g.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
$font = New-Object System.Drawing.Font 'Segoe UI', 10
$i = 0
foreach ($f in $files) {
    $x = ($i % $Across) * $Width; $y = [math]::Floor($i / $Across) * ($h + $label)
    $img = [System.Drawing.Image]::FromFile($f.FullName)
    $g.DrawImage($img, $x, $y, $Width, $h); $img.Dispose()
    $g.DrawString($f.BaseName, $font, [System.Drawing.Brushes]::White, $x + 4, $y + $h + 2)
    $i++
}
$sheet.Save($Out, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose(); $sheet.Dispose()
Write-Output "sheet: $($files.Count) shots -> $Out"
