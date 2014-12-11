# Must run using JES

# Main method converts a .jpg file to a .txt file
def main():
  # Open data and image files
  dataFile = open(pickAFile(), "w")
  image = makePicture(pickAFile())

  # Get height and width of image
  width = getWidth(image)
  height = getHeight(image)
  dataFile.write("ncols\t" + str(width) + "\n")
  dataFile.write("nrows\t" + str(height) + "\n")

  # Find RGB values for each pixel
  for y in range(0, height):
    for x in range(0, width):
      pixel = getPixel(image, x, y)
      pixelR = getRed(pixel)
      pixelG = getGreen(pixel)
      pixelB = getBlue(pixel)
      dataFile.write(str(pixelR) + " " + str(pixelG) + " " + str(pixelB) + " ")
    dataFile.write("\n")

  # Close data file
  dataFile.close()
