package me.petrolingus.iop;

import javafx.embed.swing.SwingFXUtils;
import javafx.scene.image.*;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;

public class Controller {

    public ImageView imageView;
    public ImageView imageView2;

    double[][] h12 = new double[][] {
            {0, -1, 0},
            {-1, 4, -1},
            {0, -1, 0}
    };

    double[][] h13 = new double[][] {
            {-1, -1, -1},
            {-1, 8, -1},
            {-1, -1, -1}
    };

    double[][] h14 = new double[][] {
            {1, -2, 1},
            {-2, 4, -2},
            {1, -2, 1}
    };

    double[][] h1 = new double[][] {
            {0, -1, 0},
            {-1, 5, -1},
            {0, -1, 0}
    };

    double[][] h2 = new double[][] {
            {-1, -1, -1},
            {-1, 9, -1},
            {-1, -1, -1}
    };

    double[][] h3 = new double[][] {
            {1, -2, -1},
            {-2, 5, -2},
            {1, -2, 1}
    };

    List<double[][]> masks = List.of(h1, h2, h3, h12, h13, h14);
    List<String> maskNames = List.of("h1", "h2", "h3", "h12", "h13", "h14");

    public void initialize() {

        Image image = new Image("/rose.jpeg");
        imageView.setImage(image);
        PixelReader pixelReader = image.getPixelReader();

        double[][] brightness = new double[512][512];
        for (int i = 0; i < 512; i++) {
            for (int j = 0; j < 512; j++) {
                brightness[i][j] = pixelReader.getColor(j, i).getBrightness();
            }
        }

        generateLinearSamples(brightness);
        generateFourierSamples(brightness);

        System.out.println("Done!");
    }

    private void generateLinearSamples(double[][] brightness) {
        for (int s = 0; s < 2; s++) {
            for (int m = 3; m < 6; m++) {
                for (int t = 0; t < 3; t++) {
                    double c = 0.6 + t * 0.1;
                    double[][] outline1 = linearOutline(brightness, masks.get(m));
                    double[][] smooth = smooth(brightness);
                    double[][] result = new double[512][512];
                    for (int i = 0; i < 512; i++) {
                        for (int j = 0; j < 512; j++) {
                            if (s == 0) {
                                result[i][j] = c * brightness[i][j] + (1.0 - c) * outline1[i][j];
                            } else {
                                result[i][j] = c * smooth[i][j] + (1.0 - c) * outline1[i][j];
                            }
                        }
                    }
                    if (s == 0) {
                        save("linear-" + maskNames.get(m) + "-" + c, getImageFromPixels(result));
                    } else {
                        save("smooth-linear-" + maskNames.get(m) + "-" + c, getImageFromPixels(result));
                    }
                }
            }
        }
    }

    private void generateFourierSamples(double[][] brightness) {

        for (int t = 0; t < 3; t++) {
            double c = 0.6 + t * 0.1;
            for (int rr = 50; rr < 250; rr += 50) {
                double[][] contour = fourierOutline(brightness, rr, false);
//                double[][] smooth = fourierOutline(brightness, rr, true);
                double[][] result = new double[512][512];
                for (int i = 0; i < 512; i++) {
                    for (int j = 0; j < 512; j++) {
                        result[i][j] = c * brightness[i][j] + (1.0 - c) * contour[i][j];
                    }
                }
                save("fourier-r" + rr + "-" + c, getImageFromPixels(result));
            }
        }

    }

    private void save(String name, Image img) {
        File outputFile = new File("C:\\Users\\petro\\Desktop\\samples\\" + name + ".png");
        BufferedImage bImage = SwingFXUtils.fromFXImage(img, null);
        try {
            ImageIO.write(bImage, "png", outputFile);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private double[][] fourierOutline(double[][] brightness, double rr, boolean isSmooth) {

        Gaussian gaussian = new Gaussian(1.0, 0, rr/2);

        Complex[][] spectrum = swap(ImageFourier.fft(brightness));
        for (int i = 0; i < 512; i++) {
            for (int j = 0; j < 512; j++) {
                double dx = (double) j - 256;
                double dy = (double) i - 256;
                double r = FastMath.sqrt(dx * dx + dy * dy);
                if (isSmooth) {
                    spectrum[i][j] = spectrum[i][j].multiply(gaussian.value(r));
                } else {
                    spectrum[i][j] = spectrum[i][j].multiply(1.0 - gaussian.value(r));
                }
            }
        }
        Complex[][] rev = ImageFourier.ifft(spectrum);

        double[][] result = new double[512][512];
        for (int i = 0; i < 512; i++) {
            for (int j = 0; j < 512; j++) {
                result[i][j] = rev[i][j].abs();
            }
        }

        return normalize(result);
    }

    private double[][] linearOutline(double[][] brightness, double[][] mask) {

        double[][] result = new double[512][512];
        for (int i = 1; i < 511; i++) {
            for (int j = 1; j < 511; j++) {

                double m00 = brightness[i - 1][j - 1] * mask[0][0];
                double m01 = brightness[i - 1][j] * mask[0][1];
                double m02 = brightness[i - 1][j + 1] * mask[0][2];

                double m10 = brightness[i][j - 1] * mask[1][0];
                double m11 = brightness[i][j] * mask[1][1];
                double m12 = brightness[i][j + 1] * mask[1][2];

                double m20 = brightness[i + 1][j - 1] * mask[2][0];
                double m21 = brightness[i + 1][j] * mask[2][1];
                double m22 = brightness[i + 1][j + 1] * mask[2][2];

                double sum = m00 + m01 + m02 + m10 + m11 + m12 + m20 + m21 + m22;

                result[i][j] = sum;
            }
        }

        return normalize2(result);
    }

    private double[][] smooth(double[][] brightness) {

        double[][] mask = new double[][] {
                {1, 1, 1},
                {1, 8, 1},
                {1, 1, 1}
        };

        double[][] result = new double[512][512];
        for (int i = 1; i < 511; i++) {
            for (int j = 1; j < 511; j++) {

                double m00 = brightness[i - 1][j - 1] * mask[0][0];
                double m01 = brightness[i - 1][j] * mask[0][1];
                double m02 = brightness[i - 1][j + 1] * mask[0][2];

                double m10 = brightness[i][j - 1] * mask[1][0];
                double m11 = brightness[i][j] * mask[1][1];
                double m12 = brightness[i][j + 1] * mask[1][2];

                double m20 = brightness[i + 1][j - 1] * mask[2][0];
                double m21 = brightness[i + 1][j] * mask[2][1];
                double m22 = brightness[i + 1][j + 1] * mask[2][2];

                double sum = m00 + m01 + m02 + m10 + m11 + m12 + m20 + m21 + m22;

                result[i][j] = sum;
            }
        }

        return normalize(result);
    }

    private Image getImageFromPixels(double[][] pixels) {
        int w = pixels[0].length;
        int h = pixels.length;
        int[] buffer = new int[w * h];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                int y = (int) Math.round(255 * pixels[i][j]);
                buffer[w * i + j] = 0xFF << 24 | y << 16 | y << 8 | y;
            }
        }
        WritableImage image = new WritableImage(w, h);
        PixelWriter pixelWriter = image.getPixelWriter();
        pixelWriter.setPixels(0, 0, w, h, PixelFormat.getIntArgbInstance(), buffer, 0, w);
        return image;
    }

    public double[][] normalize(double[][] matrix) {
        int w = matrix[0].length;
        int h = matrix.length;
        double max = Double.MIN_VALUE;
        double min = Double.MAX_VALUE;
        for (double[] row : matrix) {
            for (int j = 0; j < w; j++) {
                max = Math.max(max, row[j]);
                min = Math.min(min, row[j]);
            }
        }
        double[][] result = new double[h][w];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                result[i][j] = (matrix[i][j] - min) / (max - min);
            }
        }
        return result;
    }

    public double[][] normalize2(double[][] matrix) {
        int w = matrix[0].length;
        int h = matrix.length;
        double[][] result = new double[h][w];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                double value = matrix[i][j];
                if (value < 0.0) {
                    result[i][j] = 0.0;
                } else {
                    result[i][j] = Math.min(value, 1.0);
                }
            }
        }
        return result;
    }

    public static double[][] swap(double[][] matrix) {
        int w = matrix[0].length;
        int h = matrix.length;
        int hw = w / 2;
        int hh = h / 2;
        double[][] result = new double[h][w];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                int x = (j < hw) ? j + hw : j - hw;
                int y = (i < hh) ? i + h : i - hh;
                result[i][j] = matrix[y][x];
            }
        }
        return result;
    }

    public static Complex[][] swap(Complex[][] matrix) {
        int w = matrix[0].length;
        int h = matrix.length;
        int hw = w / 2;
        int hh = h / 2;
        Complex[][] result = new Complex[h][w];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                int x = (j < hw) ? j + hw : j - hw;
                int y = (i < hh) ? i + hh : i - hh;
                result[i][j] = matrix[y][x];
            }
        }
        return result;
    }
}
