// =================================================================================
// SKRIP FINAL: SENTINEL-2 BARA BULUKUMBA (ALL BANDS VERSION)
// Output: 21 Band (12 Fisik + 9 Indeks), Full Resolution 10m
// =================================================================================

// 1. INPUT DATA
var aoi = ee.FeatureCollection('projects/perairan-bara-bulukumba/assets/Peta_Lokasi_Penelitian/Daerah_Penelitian');
var lamun = ee.FeatureCollection('projects/perairan-bara-bulukumba/assets/Lamun_Interpolasi/Lamun');

// 2. FUNGSI MASKING & KOREKSI (SENTINEL-2)
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  // Scaling Sentinel-2 (Integer -> Float)
  var scaled = image.divide(10000);

  // KOREKSI NEGATIF (Clamping to 0)
  var corrected = scaled.where(scaled.lt(0), 0);

  return image.addBands(corrected, null, true)
      .updateMask(mask)
      .copyProperties(image, ["system:time_start"]);
}

// 3. FUNGSI INDEKS (9 Indeks) - NDAVI DIUBAH MENGGUNAKAN NIR & BLUE
function addS2Indices(image) {
  var eps = 0.001; 
  
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var mndwi = image.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  var ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
  
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('B8'), 'RED': image.select('B4'), 'BLUE': image.select('B2')
    }).rename('EVI');
  
  var sabi = image.expression(
    '(NIR - RED) / (BLUE + GREEN + eps)', { 
      'NIR': image.select('B8'), 'RED': image.select('B4'),
      'BLUE': image.select('B2'), 'GREEN': image.select('B3'), 'eps': eps
    }).rename('SABI');
  
  var gndvi = image.normalizedDifference(['B8', 'B3']).rename('GNDVI');
  
  // NDAVI DIUBAH: Menggunakan NIR dan Blue (B8 dan B2)
  var ndavi = image.normalizedDifference(['B8', 'B2']).rename('NDAVI');
  
  var awei = image.expression(
    '4 * (GREEN - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)', {
      'GREEN': image.select('B3'), 'SWIR1': image.select('B11'),
      'NIR': image.select('B8'), 'SWIR2': image.select('B12')
    }).rename('AWEI');
  
  var bathymetry = image.expression(
    'BLUE / (GREEN + eps)', { 
      'BLUE': image.select('B2'), 'GREEN': image.select('B3'), 'eps': eps
    }).rename('BATHYMETRY');
  
  return image.addBands([ndvi, mndwi, ndwi, evi, sabi, gndvi, ndavi, awei, bathymetry]);
}

// 4. PROSES DATA
var s2Collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate('2023-07-01', '2023-09-30') 
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(maskS2clouds)
  .map(addS2Indices);

// =================================================================================
// DEFINISI SEMUA BAND (TOTAL 21 BAND)
// =================================================================================
// 12 Band Fisik (Termasuk Red Edge yang tidak ada di Landsat)
var physicalBands = [
  'B1',  // Coastal Aerosol
  'B2',  // Blue
  'B3',  // Green
  'B4',  // Red
  'B5',  // Red Edge 1 (Bagus untuk klorofil)
  'B6',  // Red Edge 2
  'B7',  // Red Edge 3
  'B8',  // NIR
  'B8A', // Red Edge 4 (Narrow NIR)
  'B9',  // Water Vapour
  'B11', // SWIR 1
  'B12'  // SWIR 2
];

// 9 Band Indeks
var indexBands = [
  'NDVI', 'MNDWI', 'NDWI', 'EVI', 'SABI', 'GNDVI', 'NDAVI', 'AWEI', 'BATHYMETRY'
];

// Gabungkan Semua (21 Band)
var allBands = physicalBands.concat(indexBands);

var exportImage = s2Collection.median()
  .select(allBands) // Pilih SEMUA band
  .float() 
  .clip(aoi);

// =================================================================================
// 5. TAMPILAN HASIL (FULL CHECK)
// =================================================================================

// A. STATISTIK LUASAN
print('=== 1. STATISTIK LUASAN AREA ===');
var area = aoi.geometry().area();
print('Luas (Hektar):', area.divide(10000));

// B. CEK TIPE DATA LENGKAP (21 BAND)
print('=== 2. CEK TIPE DATA SIAP EXPORT ===');

print('--- A. Tipe Data 12 Band Fisik (Termasuk Red Edge) ---');
print(exportImage.select(physicalBands).bandTypes()); 

print('--- B. Tipe Data 9 Indeks ---');
print(exportImage.select(indexBands).bandTypes()); 

// C. VALIDASI KOREKSI NEGATIF
print('=== 3. VALIDASI KOREKSI NEGATIF (NILAI MINIMAL) ===');
var bandStats = exportImage.select(['B2', 'B3', 'B4', 'B8'])
  .reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: aoi.geometry(),
    scale: 10, 
    maxPixels: 1e13
  });
print('Nilai Min/Max (Min harus 0):', bandStats);

// =================================================================================
// 6. VISUALISASI & EXPORT
// =================================================================================
Map.centerObject(aoi, 14);

// 1. Layer Citra Asli (True Color)
Map.addLayer(exportImage, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3}, 'Sentinel-2 - True Color RGB');

// 2. Layer Batas Area Studi
var outline = ee.Image().byte().paint({
  featureCollection: aoi,
  color: 1,
  width: 3
});
Map.addLayer(outline, {palette: 'red'}, 'Batas Area Studi');

Export.image.toDrive({
  image: exportImage,
  description: 'Sentinel2_Bara_2023_AllBands',
  folder: 'Perairan_Bara_Bulukumba',
  fileNamePrefix: 'S2_Bara_2023_21Bands_Full',
  scale: 10, // Resolusi 10m
  region: aoi.geometry(),
  fileFormat: 'GeoTIFF',
  maxPixels: 1e13
});