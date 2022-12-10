package me.petrolingus.iop;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

public class ImageFourier {

    public static final FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);

    public static Complex[][] fft(double[][] image) {

        int w = image[0].length;
        int h = image.length;

        Complex[][] temp = new Complex[h][w];
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                temp[i][j] = new Complex(image[i][j]);
            }
        }

        return fft(temp);
    }

    public static Complex[][] fft(Complex[][] image) {
        return transform(image, TransformType.FORWARD);
    }

    public static Complex[][] ifft(Complex[][] image) {
        return transform(image, TransformType.INVERSE);
    }

    private static Complex[][] transform(Complex[][] image, TransformType transformType) {

        int w = image[0].length;
        int h = image.length;

        Complex[][] temp1 = new Complex[h][w];
        for (int i = 0; i < h; i++) {
            temp1[i] = fft.transform(image[i], transformType);
        }

        Complex[][] temp2 = new Complex[h][w];
        for (int i = 0; i < h; i++) {
            Complex[] col = new Complex[h];
            for (int j = 0; j < h; j++) {
                col[j] = temp1[j][i];
            }
            Complex[] transform = fft.transform(col, transformType);
            for (int j = 0; j < h; j++) {
                temp2[j][i] = transform[j];
            }
        }

        return temp2;
    }
}
