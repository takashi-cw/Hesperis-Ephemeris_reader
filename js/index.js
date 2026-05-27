/**
 * stella-bsp-reader — JS エントリポイント
 *
 * JPL DE天体暦（.bsp / SPK Type 2）を Node.js またはブラウザから
 * ビルドツールなしで使える純 JavaScript ライブラリ。
 */

export { loadBsp, parseBsp, BspFile } from './src/bsp-reader.js';
export { getCoverageJd, formatCoverageMessage, assertInCoverage } from './src/bsp-validator.js';
export { NAIF, PLANET_NAIF, J2000_JD, AU_KM } from './src/constants.js';
