# PowerShell version of the version generation script
param(
    [string]$OutputPath = ""
)

$ErrorActionPreference = "Stop"

function Get-GitRevisionHash {
    try {
        $hash = git rev-parse --short HEAD 2>$null
        return $hash
    }
    catch {
        return "none"
    }
}

function Gen-VersionGit {
    $versionTemplate = @"
 character(LEN=32),parameter :: schism_version = '@{VERSION_SCHISM}', git_version = '@{VERSION_GIT}'
"@

    $versionPath = Split-Path -Parent $PSCommandPath
    $versionScratchPath = Join-Path $versionPath "_version"
    $describeRegex = '^(v?[0-9]+\.[0-9]+\.[0-9]+)(-[0-9]+)?(-g[a-z0-9]+)?$'

    try {
        git describe --always --dirty | Out-File -FilePath $versionScratchPath -Encoding ASCII -NoNewline 2>$null
    }
    catch {
        "unknown" | Out-File -FilePath $versionScratchPath -Encoding ASCII -NoNewline
    }

    $versionRaw = (Get-Content $versionScratchPath -Raw).Trim()
    Remove-Item $versionScratchPath -Force -ErrorAction SilentlyContinue

    $isDirty = $false
    if ($versionRaw -like "*-dirty") {
        $isDirty = $true
        $versionRaw = $versionRaw -replace "-dirty$", ""
    }

    if ($versionRaw -match $describeRegex) {
        $schismVersion = $Matches[1]
        $gitVersion = $Matches[3] -replace "^-g", ""
        
        if ($Matches[2] -or $isDirty) {
            $nmod = $Matches[2] -replace "^-", ""
            $schismVersion += "mod"
            $gitVersion += " ($nmod commits since semantic tag, edits=$isDirty)"
        }
    }
    else {
        $schismVersion = "semantic version not determined"
        $gitVersion = Get-GitRevisionHash
    }

    return @($gitVersion, $schismVersion)
}

function Gen-VersionUser {
    $schismUserVersionFile = "schism_version_user.txt"
    $userVersionPath = Join-Path (Split-Path -Parent $PSCommandPath) $schismUserVersionFile
    
    if (Test-Path $userVersionPath) {
        $userVersion = (Get-Content $userVersionPath -Raw).Trim()
        if ($userVersion.Length -lt 3) {
            $userVersion = "develop"
        }
    }
    else {
        $userVersion = "develop"
    }
    
    return $userVersion
}

function Gen-Version {
    param([string]$OutputFile = "")
    
    $versionPath = Split-Path -Parent $PSCommandPath
    $templatePath = Join-Path $versionPath "schism_version.F90.template"
    
    if (-not $OutputFile) {
        $versionFilePath = Join-Path $versionPath "schism_version.F90"
    }
    else {
        $versionFilePath = $OutputFile
    }
    
    $gitVersion, $schismVersion = Gen-VersionGit
    
    if ($schismVersion -eq "semantic version not determined") {
        $schismVersion = Gen-VersionUser
    }
    
    if (-not $gitVersion -or $gitVersion -eq "none") {
        $gitVersion = Get-GitRevisionHash
    }
    
    Write-Host " SCHISM version: $schismVersion"
    Write-Host " GIT commit $gitVersion"
    
    $templateContent = Get-Content $templatePath -Raw
    $outputContent = $templateContent -replace '@\{VERSION_GIT\}', $gitVersion -replace '@\{VERSION_SCHISM\}', $schismVersion
    
    $outputContent | Out-File -FilePath $versionFilePath -Encoding ASCII
}

# Main execution
if ($MyInvocation.InvocationName -ne '.') {
    Gen-Version $OutputPath
}