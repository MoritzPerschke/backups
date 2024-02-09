#include "image_diff.h"
/* 
Compute PSNR (https://www.ni.com/de-at/shop/data-acquisition-and-control/add-ons-for-data-acquisition-and-control/what-is-vision-development-module/peak-signal-to-noise-ratio-as-an-image-quality-metric.html) 
*/
double computePSNR(const Image& original, const Image& compare){
    int numPixels = original.height * original.width;
    double mse = 0.;

    for (int i = 0; i < numPixels; i++){
        Color orig = original.pixels[i];
        Color comp = compare.pixels[i];
        double distSquare = pow((orig.x - comp.x), 2)
                          + pow((orig.y - comp.y), 2)
                          + pow((orig.z - comp.z), 2);
        mse += distSquare;
    }
    mse /= numPixels;

    return 20 * log10(pow(255.0, 2) / mse);
}
int main(int argc, char* argv[]){

    if(argc != 3){
        std::cerr << err << "Usage: <samples of image> <sampler of image>";
        return 1;
    }

    int samples = atoi(argv[1]);
    Rand_Gen_t gen = Rand_Gen_t(atoi(argv[2]));
    
    // This path is used to run the program from 'Solution/build/bin'
    const char* original_path = "images/original.ppm";
    const std::string image2_path = "images/" + to_string(gen) + "/" + std::to_string(samples) + ".ppm";
    // std::cout << ok << "image paths are: " << original_path << " and " << image2_path << std::endl;

    Image original = Image(original_path);
    Image compare = Image(image2_path);

    if (original.width != compare.width || original.height != compare.height){
        std::cout << err << "Image dimensions do not match, exiting..." << std::endl;
        return 1;
    }

    double psnr = computePSNR(original, compare);
    std::cerr << "{\"psnr\":" << psnr << "}" << std::endl;

    /* Calculate difference image */
    Image diff(original.width, original.height);
    for (int y = 0; y < original.height; y++){
        for (int x = 0; x < original.width; x++){
            Color a = original.getColor(x, y);
            Color b = compare.getColor(x, y);
            diff.setColor(x, y, sub_abs(a, b));
        }
    }
    std::cout << ok << "Calculation Done, saving..." << std::endl;
    diff.SaveRaw("images/" + to_string(gen) +  "/differences/diff_" + std::to_string(samples) + ".ppm");

    return 0;
}