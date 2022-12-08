package me.petrolingus.iop;

import javafx.embed.swing.SwingFXUtils;
import javafx.scene.image.*;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.util.FastMath;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;

public class Controller {

    public ImageView imageView;
    public ImageView imageView2;

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

        double c = 0.5;

        double[][] outline1 = outline2(brightness);

        double[][] result = new double[512][512];
        for (int i = 0; i < 512; i++) {
            for (int j = 0; j < 512; j++) {
                result[i][j] = c * brightness[i][j] + (1.0 - c) * outline1[i][j];
            }
        }

        imageView2.setImage(getImageFromPixels(result));

        save("234234432", getImageFromPixels(result));
    }

    private void save(String name, Image img) {

//        File outputFile = new File("C:\\samples\\" + name + ".jpg");
        File outputFile = new File(name + ".jpg");
        BufferedImage bImage = SwingFXUtils.fromFXImage(img, null);
        try {
            ImageIO.write(bImage, "jpg", outputFile);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }


    }

    private double[][] outline1(double[][] brightness) {

        Complex[][] spectrum = swap(ImageFourier.fft(brightness));
        for (int i = 0; i < 512; i++) {
            for (int j = 0; j < 512; j++) {
                double dx = (double) j - 256;
                double dy = (double) i - 256;
                double r = FastMath.sqrt(dx * dx + dy * dy);
                if (r < 10) {
                    spectrum[i][j] = Complex.ZERO;
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

    private double[][] outline2(double[][] brightness) {

        double[][] mask = new double[][] {
                {-1, -1, -1},
                {-1, 8, -1},
                {-1, -1, -1}
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

                if (sum < 0) {
                    sum = 0;
                } else if (sum > 1) {
                    sum = 1;
                }

                result[i][j] = sum;
            }
        }

        return result;
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
