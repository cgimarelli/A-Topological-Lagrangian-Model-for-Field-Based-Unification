import PyInstaller.__main__

print("Starting build... please wait...")

PyInstaller.__main__.run([
    'SpaceScanner.py',
    '--name=SpaceScanner',
    '--onefile',
    '--clean',
    '--noconfirm',
    '--noconsole',
    '--icon=SpaceScannerIcon.ico',
    '--collect-all=customtkinter',
    '--collect-all=astroquery',    
])