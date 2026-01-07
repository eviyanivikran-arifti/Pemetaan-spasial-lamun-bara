// =================================================================================
// SKRIP FINAL: LANDSAT 8 BARA BULUKUMBA
// Output: 16 Band, Cek Tipe Data Lengkap, Visualisasi: True Color & Area Studi
// =================================================================================

// 1. INPUT DATA
var aoi = ee.FeatureCollection('projects/perairan-bara-bulukumba/assets/Peta_Lokasi_Penelitian/Daerah_Penelitian');
var lamun = ee.FeatureCollection('projects/perairan-bara-bulukumba/assets/Lamun_Interpolasi/Lamun');

// 2. FUNGSI MASKING & KOREKSI
function maskL8clouds(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0).and(qa.bitwiseAnd(1 << 4).eq(0));
  
  // Scaling
  var scaled = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  
  // KOREKSI NEGATIF (Clamping to 0) -> Agar Min value jadi 0
  var corrected = scaled.where(scaled.lt(0), 0);
  
  return image.addBands(corrected, null, true)
      .updateMask(mask)
      .copyProperties(image, ["system:time_start"]);
}

// 3. FUNGSI INDEKS (9 Indeks) - NDAVI DIUBAH MENGGUNAKAN NIR & BLUE
function addL8Indices(image) {
  var eps = 0.001; 
  var ndvi = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  var mndwi = image.normalizedDifference(['SR_B3', 'SR_B6']).rename('MNDWI');
  var ndwi = image.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('SR_B5'), 'RED': image.select('SR_B4'), 'BLUE': image.select('SR_B2')
    }).rename('EVI');
  var sabi = image.expression(
    '(NIR - RED) / (BLUE + GREEN + eps)', { 
      'NIR': image.select('SR_B5'), 'RED': image.select('SR_B4'),
      'BLUE': image.select('SR_B2'), 'GREEN': image.select('SR_B3'), 'eps': eps
    }).rename('SABI');
  var gndvi = image.normalizedDifference(['SR_B5', 'SR_B3']).rename('GNDVI');
  
  // NDAVI DIUBAH: Menggunakan NIR dan Blue (SR_B5 dan SR_B2)
  var ndavi = image.normalizedDifference(['SR_B5', 'SR_B2']).rename('NDAVI');
  
  var awei = image.expression(
    '4 * (GREEN - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)', {
      'GREEN': image.select('SR_B3'), 'SWIR1': image.select('SR_B6'),
      'NIR': image.select('SR_B5'), 'SWIR2': image.select('SR_B7')
    }).rename('AWEI');
  var bathymetry = image.expression(
    'BLUE / (GREEN + eps)', { 
      'BLUE': image.select('SR_B2'), 'GREEN': image.select('SR_B3'), 'eps': eps
    }).rename('BATHYMETRY');
  
  return image.addBands([ndvi, mndwi, ndwi, evi, sabi, gndvi, ndavi, awei, bathymetry]);
}

// 4. PROSES DATA
var l8Collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(aoi)
  .filterDate('2023-07-01', '2023-09-30')
  .filter(ee.Filter.lt('CLOUD_COVER', 20))
  .map(maskL8clouds)
  .map(addL8Indices);

// Pilih Band untuk Export (16 Band Total)
var bandNames = [
  'SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', // 7 Fisik
  'NDVI', 'MNDWI', 'NDWI', 'EVI', 'SABI', 'GNDVI', 'NDAVI', 'AWEI', 'BATHYMETRY' // 9 Indeks
];

var exportImage = l8Collection.median()
  .select(bandNames)
  .float() // Paksa tipe data Float
  .clip(aoi);

// =================================================================================
// 5. TAMPILAN HASIL (FULL CHECK)
// =================================================================================

// A. STATISTIK LUASAN
print('=== 1. STATISTIK LUASAN AREA ===');
var area = aoi.geometry().area();
print('Luas (Hektar):', area.divide(10000));

// B. CEK TIPE DATA LENGKAP
print('=== 2. CEK TIPE DATA SIAP EXPORT ===');

// Tampilkan Band Fisik
print('--- A. Tipe Data 7 Band Fisik (Harus "float") ---');
var physicalTypes = exportImage.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']).bandTypes();
print(physicalTypes); 

// Tampilkan Indeks
print('--- B. Tipe Data 9 Indeks (Harus "float") ---');
var allIndexNames = ['NDVI', 'MNDWI', 'NDWI', 'EVI', 'SABI', 'GNDVI', 'NDAVI', 'AWEI', 'BATHYMETRY'];
var indexTypes = exportImage.select(allIndexNames).bandTypes();
print(indexTypes); 

// C. VALIDASI KOREKSI NEGATIF
print('=== 3. VALIDASI KOREKSI NEGATIF (NILAI MINIMAL) ===');
var bandStats = exportImage.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5'])
  .reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: aoi.geometry(),
    scale: 30,
    maxPixels: 1e13
  });
print('Nilai Min/Max (Min harus 0):', bandStats);

// =================================================================================
// 6. VISUALISASI & EXPORT
// =================================================================================
Map.centerObject(aoi, 14);

// 1. Layer Citra Asli (True Color)
Map.addLayer(exportImage, {bands: ['SR_B4', 'SR_B3', 'SR_B2'], min: 0, max: 0.3}, 'Landsat 8 - True Color RGB');

// 2. Layer Batas Area Studi (Garis Merah, Tengah Kosong)
var empty = ee.Image().byte();
var outline = empty.paint({
  featureCollection: aoi,
  color: 1,
  width: 3
});
Map.addLayer(outline, {palette: 'red'}, 'Batas Area Studi');

Export.image.toDrive({
  image: exportImage,
  description: 'Landsat8_Bara_2023_Ready',
  folder: 'Perairan_Bara_Bulukumba',
  fileNamePrefix: 'L8_Bara_2023_16Bands_Final',
  scale: 30,
  region: aoi.geometry(),
  fileFormat: 'GeoTIFF',
  maxPixels: 1e13
});