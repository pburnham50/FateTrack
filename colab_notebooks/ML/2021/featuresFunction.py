def greyprops_far(regionmask, intensity):
    return greycoprops(greycomatrix(intensity[regionmask], [10], angles=[0], levels=256, symmetric=False, normed=True),
                       prop=('contrast', 'dissimilarity', 'correlation', 'energy', 'homogeneity')

def greyprops_mid(regionmask, intensity):
    return greycoprops(greycomatrix(intensity[regionmask], [5], angles=[0], levels=256, symmetric=False, normed=True),
                       prop=('contrast', 'dissimilarity', 'correlation', 'energy', 'homogeneity')

def greyprops_near(regionmask, intensity):
    return greycoprops(greycomatrix(intensity[regionmask], distances = [2], angles=[0], levels=256, symmetric=False, normed=True),
                       prop=('contrast', 'dissimilarity', 'correlation', 'energy', 'homogeneity')

def spots(regionmask, intensity):
    spotlabels = measure.label(intensity[regionmask], connectivity=2, background=0)
    npspots = np.mask(np.unique(spotlabels))
    spotpixels = len(spotlabels[spotlabels!=0])
    return (npspots, spotpixels)

def imvar(regionmask, intensity):
    return np.sum((intensity[regionmask]-np.mean(intensity[regionmask]))**2)
                       
def nuclearFeatureTable(PathToImage, subImageIndex, round):
  segFile = PathToImage + 'masks/sub_' + str(subImageIndex) + '_' + str(round) + 'DAPI_seg.npy'
  newnp = np.asarray(np.load(segFile, allow_pickle=True)).item()
  outs = (newnp['outlines']>0).astype(int)
  masks = newnp['masks']
  imorig = newnp['img']

  image_props = measure.regionprops_table(masks, intensity_image=imorig,
                                          properties=('label','area','filled_area', 'bbox', 'centroid',
                                                      'eccentricity','solidity','convex_area',
                                                      'mean_intensity','min_intensity','max_intensity',
                                                      'orientation','major_axis_length','minor_axis_length',
                                                      'perimeter','extent','intensity_image'),
                                          extra_properties = (greyprops_near, greyprops_mid, greyprops_far, spots, imvar))
  im_df = pd.DataFrame(image_props)
  export = im_df.to_csv(PathToImage + 'tables/sub_' + str(subImageIndex) + "_" + str(round) + "_nuclei.feature.csv", index = None, header=True)


                       
