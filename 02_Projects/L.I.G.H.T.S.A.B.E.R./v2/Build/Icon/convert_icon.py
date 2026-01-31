from PIL import Image

# Open the PNG
img = Image.open("Space_Scanner_Icon.png")

# Save as ICO with multiple sizes (windows likes 256x256, 128x128, etc)
img.save("SpaceScannerIcon.ico", format='ICO', sizes=[(256, 256), (128, 128), (64, 64), (48, 48), (32, 32), (16, 16)])

print("Converted icon.png to SpaceScanner.ico!")