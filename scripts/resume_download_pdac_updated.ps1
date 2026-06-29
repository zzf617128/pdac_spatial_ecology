$ErrorActionPreference = "Stop"

$Url = "https://zenodo.org/api/records/10712047/files/PDAC_Updated.rds/content"
$OutDir = "E:\PDAC_TLS\pdac_spatial_ecology\data\external\GSE272362_zenodo"
$PartPath = Join-Path $OutDir "PDAC_Updated.rds.part"
$FinalPath = Join-Path $OutDir "PDAC_Updated.rds"
$ExpectedBytes = 11638143468

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

if (Test-Path -LiteralPath $FinalPath) {
    $existing = Get-Item -LiteralPath $FinalPath
    if ($existing.Length -eq $ExpectedBytes) {
        Write-Host "Already complete: $FinalPath"
        exit 0
    }
    Write-Host "Existing final file has unexpected size; resuming into .part instead."
}

Write-Host "Downloading with resume support..."
Write-Host "URL: $Url"
Write-Host "Output: $PartPath"

curl.exe `
  -L `
  --fail `
  --retry 30 `
  --retry-delay 10 `
  --retry-all-errors `
  --connect-timeout 60 `
  --continue-at - `
  --output "$PartPath" `
  "$Url"

$downloaded = Get-Item -LiteralPath $PartPath
Write-Host "Downloaded bytes: $($downloaded.Length) / $ExpectedBytes"

if ($downloaded.Length -eq $ExpectedBytes) {
    Move-Item -Force -LiteralPath $PartPath -Destination $FinalPath
    Write-Host "Download complete: $FinalPath"
} else {
    Write-Host "Download incomplete. Re-run this script later to resume."
    exit 2
}

