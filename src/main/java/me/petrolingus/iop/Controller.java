package me.petrolingus.iop;

import javafx.scene.control.Slider;
import javafx.scene.control.TextField;
import javafx.scene.image.*;

import java.util.List;

public class Controller {

    public ImageView imageView;
    public ImageView imageView2;

    public TextField h00;
    public TextField h01;
    public TextField h02;
    public TextField h10;
    public TextField h11;
    public TextField h12;
    public TextField h20;
    public TextField h21;
    public TextField h22;

    public Slider radiusSlider;
    public Slider mixingSlider;

    double[][] brightness;
    double[][] result;
    String method = null;

    public void initialize() {

        Image image = new Image("/forest.jpg");

        int width = (int) image.getWidth();
        int height = (int) image.getHeight();
        brightness = new double[width][height];

        PixelReader pixelReader = image.getPixelReader();
        for (int i = 0; i < 512; i++) {
            for (int j = 0; j < 512; j++) {
                brightness[i][j] = pixelReader.getColor(j, i).getBrightness();
            }
        }

        imageView.setImage(getImageFromPixels(brightness));

        FooBar fooBar = new FooBar();
        fooBar.generateLinearSamples(brightness);
        fooBar.generateFourierSamples(brightness);
    }

    public void onTimeProcessing() {

        method = "linear";

        double h00 = Double.parseDouble(this.h00.getText());
        double h01 = Double.parseDouble(this.h01.getText());
        double h02 = Double.parseDouble(this.h02.getText());
        double h10 = Double.parseDouble(this.h10.getText());
        double h11 = Double.parseDouble(this.h11.getText());
        double h12 = Double.parseDouble(this.h12.getText());
        double h20 = Double.parseDouble(this.h20.getText());
        double h21 = Double.parseDouble(this.h21.getText());
        double h22 = Double.parseDouble(this.h22.getText());

        double[][] mask = new double[][] {
                {h00, h01, h02},
                {h10, h11, h12},
                {h20, h21, h22}
        };

        double[][] processed = new double[brightness.length][brightness.length];
        for (int i = 1; i < brightness.length - 1; i++) {
            for (int j = 1; j < brightness.length - 1; j++) {
                double m00 = brightness[i - 1][j - 1] * mask[0][0];
                double m01 = brightness[i - 1][j] * mask[0][1];
                double m02 = brightness[i - 1][j + 1] * mask[0][2];
                double m10 = brightness[i][j - 1] * mask[1][0];
                double m11 = brightness[i][j] * mask[1][1];
                double m12 = brightness[i][j + 1] * mask[1][2];
                double m20 = brightness[i + 1][j - 1] * mask[2][0];
                double m21 = brightness[i + 1][j] * mask[2][1];
                double m22 = brightness[i + 1][j + 1] * mask[2][2];
                processed[i][j] = m00 + m01 + m02 + m10 + m11 + m12 + m20 + m21 + m22;
            }
        }

        double mixingRatio = mixingSlider.getValue();

        result = new double[brightness.length][brightness.length];
        for (int i = 0; i < brightness.length; i++) {
            for (int j = 0; j < brightness.length; j++) {
                result[i][j] = mixingRatio * brightness[i][j] + (1.0 - mixingRatio) * processed[i][j];
            }
        }

        imageView2.setImage(getImageFromPixels(normalize2(result)));
    }

    public void onFrequencyProcessing() {

    }

    public void onSave() {

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

}
