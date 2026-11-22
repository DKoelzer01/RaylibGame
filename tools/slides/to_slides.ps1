# Shots into menu slides (#1071): every PNG in -Dir scaled to -Width px wide (1600) and saved as JPEG (quality 88) into
# -Out under the same name. (The game's --to-slide does one file the same way but opens the game; this needs no slot.)
#   powershell -File tools/slides/to_slides.ps1 -Dir <shots> -Out resources/menu/slides [-Width 1600] [-Quality 88]
param([string]$Dir, [string]$Out, [int]$Width = 1600, [int]$Quality = 88)
Add-Type -AssemblyName System.Drawing
$codec = [System.Drawing.Imaging.ImageCodecInfo]::GetImageEncoders() | Where-Object { $_.MimeType -eq 'image/jpeg' }
$params = New-Object System.Drawing.Imaging.EncoderParameters 1
$params.Param[0] = New-Object System.Drawing.Imaging.EncoderParameter ([System.Drawing.Imaging.Encoder]::Quality), ([long]$Quality)
New-Item -ItemType Directory -Force $Out | Out-Null
foreach ($f in Get-ChildItem -Path $Dir -Filter *.png | Sort-Object Name) {
    $img = [System.Drawing.Image]::FromFile($f.FullName)
    $w = [math]::Min($Width, $img.Width); $h = [int]($img.Height * $w / $img.Width)
    $bmp = New-Object System.Drawing.Bitmap $w, $h
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $g.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::HighQualityBicubic
    $g.DrawImage($img, 0, 0, $w, $h)
    $dest = Join-Path $Out ($f.BaseName + '.jpg')
    $bmp.Save($dest, $codec, $params)
    $g.Dispose(); $bmp.Dispose(); $img.Dispose()
    Write-Output ("slide: {0} ({1}x{2}, {3:N0} KB)" -f $dest, $w, $h, ((Get-Item $dest).Length / 1KB))
}
