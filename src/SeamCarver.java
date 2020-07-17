import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.StdOut;
import edu.princeton.cs.algs4.Stopwatch;

public class SeamCarver {
    private static final int ENERGY_MAX = 1000;
    private Picture picture;
    private double[][] energyGrid;
    private int[][] fromPixel;


    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) {
            throw new IllegalArgumentException();
        }
        this.picture = new Picture(picture);
        energyGrid = new double[height()][width()];
        fromPixel = new int[height()][width()];

        for (int k = 0; k < height(); k++) {
            for (int m = 0; m < width(); m++) {
                energyGrid[k][m] = energy(m, k);
            }
        }
    }

    // current picture
    public Picture picture() {
        return new Picture(picture);
    }

    // width of current picture
    public int width() {
        return picture.width();
    }

    // height of current picture
    public int height() {
        return picture.height();
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || y < 0 || x > width() - 1 || y > height() - 1) {
            throw new IllegalArgumentException();
        }
        if (x == 0 || y == 0 || x == width() - 1 || y == height() - 1) {
            return ENERGY_MAX;
        } else {
            return Math.sqrt(xEnergy(x, y) + yEnergy(x, y));
        }
    }

    private double xEnergy(int x, int y) {
        int leftRGB = picture.getRGB(x - 1, y);

        int leftRed = leftRGB >> 16 & 0xFF;
        int leftGreen = leftRGB >> 8 & 0xFF;
        int leftBlue = leftRGB & 0xFF;

        int rightRGB = picture.getRGB(x + 1, y);
        int rightRed = rightRGB >> 16 & 0xFF;
        int rightGreen = rightRGB >> 8 & 0xFF;
        int rightBlue = rightRGB & 0xFF;

        return (leftRed - rightRed) * (leftRed - rightRed)
                + (leftGreen - rightGreen) * (leftGreen - rightGreen) +
                (leftBlue - rightBlue) * (leftBlue - rightBlue);
    }

    private double yEnergy(int x, int y) {
        int bottomRGB = picture.getRGB(x, y + 1);
        int bottomRed = bottomRGB >> 16 & 0xFF;
        int bottomGreen = bottomRGB >> 8 & 0xFF;
        int bottomBlue = bottomRGB & 0xFF;

        int topRGB = picture.getRGB(x, y - 1);
        int topRed = topRGB >> 16 & 0xFF;
        int topGreen = topRGB >> 8 & 0xFF;
        int topBlue = topRGB & 0xFF;

        return (topRed - bottomRed) * (topRed - bottomRed)
                + (topGreen - bottomGreen) * (topGreen - bottomGreen) +
                (topBlue - bottomBlue) * (topBlue - bottomBlue);
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        transpose();
        int[] horizontalSeam = findVerticalSeam();
        transpose();
        return horizontalSeam.clone();

    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        if (height() == 1) {
            return new int[]{width() - 1};
        }
        double[][] pathTo = new double[height()][width()]; // the current minimum distance to each pixel in the picture
        for (int row = 1; row < height(); row++) {
            for (int col = 0; col < width(); col++) {
                findVertical(row, col, pathTo); // find the minimum distances to each pixel
            }
        }
        // start from the bottom row and travel up the seam
        double minEnergyLastRow = Double.POSITIVE_INFINITY;
        int lastCol = 0; // the column of the pixel with lowest pathTo value
        for (int i = 0; i < width(); i++) {
            if (pathTo[height() - 1][i] < minEnergyLastRow) { // find the lowest energy pixel in the last row
                minEnergyLastRow = pathTo[height() - 1][i];
                lastCol = i;
            }
        }
        int[] verticalSeam = new int[height()];
        for (int i = height() - 1; i > 0; i--) { // start from bottom row and go up
            verticalSeam[i] = lastCol; // each array entry corresponds to the seam column value in row [i];
            lastCol = fromPixel[i][lastCol]; // update the lastCol entry to the parent pixel
        }
        verticalSeam[0] = verticalSeam[1];
        return verticalSeam.clone();
    }


    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        transpose();
        validateSeam(seam);
        removeVerticalSeam(seam);
        transpose();
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        validateSeam(seam);

        Picture verticallyCarved = new Picture(width() - 1, height());
        for (int i = 0; i < height(); i++) {
            for (int j = 0; j < width() - 1; j++) {
                if (j < seam[i]) {
                    verticallyCarved.set(j, i, picture.get(j, i));
                } else {
                    verticallyCarved.set(j, i, picture.get(j + 1, i));
                }
            }
        }
        picture = verticallyCarved;
        energyGrid = new double[height()][width()];
        for (int k = 0; k < height(); k++) {
            for (int m = 0; m < width(); m++) {
                energyGrid[k][m] = energy(m, k);
            }
        }
    }

    // helper method to update the fromPixel[][] field by
    // calculating the minimum distance to input pixel
    private void findVertical(int row, int col, double[][] pathTo) {
        int minCol = col; // the column of the lowest-distance pixel in the row above the input pixel;
        double currentMin = Double.POSITIVE_INFINITY; // current value of nearest distance to pixel in possible paths above

        if (row == 0) {
            pathTo[row][col] = energyGrid[row][col]; // top row distances are automatically initialized to their energy value
            // currentMin = pathTo[row][col];
        } else { // initialize every other pixel to distance infinity
            pathTo[row][col] = Double.POSITIVE_INFINITY;
        }

        if (row == 1) { // second row values are always 10000 + current energy value
            pathTo[row][col] = energyGrid[row][col] + energyGrid[row - 1][col];
            currentMin = pathTo[row - 1][col];
        }
        // default path to pixel is directly above
        if (row > 1) {
            currentMin = pathTo[row - 1][col];
            minCol = col;
        }

        // check if upper left pixel is lower distance
        if (col > 1 && row > 1) {
            if (pathTo[row - 1][col - 1] < currentMin) {
                currentMin = pathTo[row - 1][col - 1];
                minCol = col - 1;
            }
        }

        // check if upper right pixel is lower distance
        if (col != width() - 1 && row > 1) {
            if (pathTo[row - 1][col + 1] < currentMin) {
                currentMin = pathTo[row - 1][col + 1];
                minCol = col + 1;
            }
        }

        // update the new path
        if (energyGrid[row][col] + currentMin < pathTo[row][col]) {
            pathTo[row][col] = energyGrid[row][col] + currentMin;
            fromPixel[row][col] = minCol; // holds the column of the nearest pixel in the above row
        }
    }

    private void transpose() {
        Picture transposed = new Picture(height(), width()); // swap the dimensions of the original picture
        double[][] newEnergy = new double[width()][height()];
        for (int i = 0; i < width(); i++) {
            for (int j = 0; j < height(); j++) {
                transposed.set(j, i, picture.get(i, j));
                newEnergy[i][j] = energyGrid[j][i];
            }
        }
        energyGrid = newEnergy;
        picture = transposed;
        fromPixel = new int[height()][width()];
    }

    private void validateSeam(int[] seam) {
        if (seam == null) {
            throw new IllegalArgumentException();
        }
        if (width() <= 1) {
            throw new IllegalArgumentException();
        }
        if (seam.length != picture.height()) {
            throw new IllegalArgumentException();
        }
        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] > picture.width() - 1) {
                throw new IllegalArgumentException();
            }
        }
        for (int i = 1; i < seam.length; i++) {
            if (Math.abs(seam[i] - seam[i - 1]) > 1) {
                throw new IllegalArgumentException();
            }
        }

    }

    //  unit testing (optional)
    public static void main(String[] args) {
// Command line input of image to be used and cols + rows to be removed
        SeamCarver picCarve = new SeamCarver(new Picture(args[0]));
        int removeCols = Integer.parseInt(args[1]);
        int removeRows = Integer.parseInt(args[2]);

        // Show original image, width, height, and energy of top left pixel
        picCarve.picture().show();
        StdOut.println("Width: " + picCarve.width());
        StdOut.println("Height: " + picCarve.height());
        StdOut.println("Energy at (0,0): " + picCarve.energy(0, 0));

        // Remove rows and begin timing
        StdOut.println("Removing " + removeCols + " columns");
        StdOut.println("Removing " + removeRows + " rows");
        Stopwatch sw = new Stopwatch();

        // Begin to find and remove specified number of verti. and horiz. seams
        for (int i = 0; i < removeCols; i++) {
            int[] columnSeam = picCarve.findVerticalSeam();
            picCarve.removeVerticalSeam(columnSeam);
        }

        for (int j = 0; j < removeRows; j++) {
            int[] rowSeam = picCarve.findHorizontalSeam();
            picCarve.removeHorizontalSeam(rowSeam);
        }

        // Show seam-carved image
        StdOut.println("Resizing time: " + sw.elapsedTime() + " seconds.");
        picCarve.picture().show();
    }

}
