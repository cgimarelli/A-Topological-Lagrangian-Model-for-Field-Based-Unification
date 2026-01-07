
import React, { useState, useCallback, useEffect } from 'react';
import ScannerForm from './components/ScannerForm';
import SkyChart from './components/SkyChart';
import AnalysisPanel from './components/AnalysisPanel';
import { runN12Scanner } from './services/scannerLogic';
import { analyzeScanResults } from './services/geminiService';
import { ScanResult, ScannerParams, AIAnalysis } from './types';

const App: React.FC = () => {
  const [params, setParams] = useState<ScannerParams>({
    start_xyz: 134000000,
    search_radius: 500000,
    base_ra: 3.542,
    base_dec: -27.80
  });

  const [results, setResults] = useState<ScanResult[]>([]);
  const [isScanning, setIsScanning] = useState(false);
  const [analysis, setAnalysis] = useState<AIAnalysis | null>(null);
  const [isAnalyzing, setIsAnalyzing] = useState(false);

  const handleScan = useCallback(async () => {
    setIsScanning(true);
    setAnalysis(null);
    
    // Simulated intergalactic hardware synchronization delay
    await new Promise(resolve => setTimeout(resolve, 1200));
    
    // Generate data
    const scanData = runN12Scanner(params);
    setResults(scanData);
    setIsScanning(false);

    if (scanData.length > 0) {
      setIsAnalyzing(true);
      try {
        const aiResponse = await analyzeScanResults(scanData);
        setAnalysis(aiResponse);
      } catch (err) {
        console.error("AI Analysis failed", err);
      } finally {
        setIsAnalyzing(false);
      }
    }
  }, [params]);

  useEffect(() => {
    handleScan();
  }, []);

  return (
    <div className={`max-w-7xl mx-auto px-4 py-8 md:py-12 transition-all duration-700 ${isScanning ? 'is-scanning pointer-events-none' : ''}`}>
      <div className="scanning-line"></div>
      
      <header className="mb-12 border-b border-slate-800 pb-8 flex flex-col md:flex-row md:items-end justify-between gap-6">
        <div className="space-y-2">
          <div className="flex items-center gap-3">
            <div className="w-10 h-10 bg-cyan-500/20 border border-cyan-500 rounded flex items-center justify-center animate-pulse">
              <svg xmlns="http://www.w3.org/2000/svg" className="h-6 w-6 text-cyan-400" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 6V4m0 2a2 2 0 100 4m0-4a2 2 0 110 4m-6 8a2 2 0 100-4m0 4a2 2 0 110-4m0 4v2m0-6V4m6 6v10m6-2a2 2 0 100-4m0 4a2 2 0 110-4m0 4v2m0-6V4" />
              </svg>
            </div>
            <h1 className="text-4xl md:text-5xl font-extrabold tracking-tighter text-white glow-text uppercase">
              <span className="text-cyan-500">DARK SECTOR</span> ANALYSIS
            </h1>
          </div>
          <p className="text-slate-400 text-xs tracking-widest uppercase opacity-70 border-l border-cyan-900 pl-4">
            Inversion Tracking &bull; Metric Tension Ensembles &bull; N12.00 Calibration Active
          </p>
        </div>
        
        <div className="flex gap-4">
          <div className="bg-slate-900/80 border border-slate-700 px-5 py-3 rounded-lg backdrop-blur-md">
            <span className="block text-[10px] text-slate-500 uppercase tracking-widest mb-1">Grid Points Analyzed</span>
            <span className="text-2xl font-bold text-cyan-400 tabular-nums">{results.length}</span>
          </div>
          <div className="bg-slate-900/80 border border-slate-700 px-5 py-3 rounded-lg backdrop-blur-md">
            <span className="block text-[10px] text-slate-500 uppercase tracking-widest mb-1">Survey Mode</span>
            <span className={`text-2xl font-bold ${params.start_xyz > 1000000 ? 'text-indigo-400' : 'text-green-400'}`}>
              {params.start_xyz > 1000000 ? 'DEEP FIELD' : 'LOCAL SKY'}
            </span>
          </div>
        </div>
      </header>

      <main className="space-y-10">
        <ScannerForm 
          params={params} 
          onChange={setParams} 
          onScan={handleScan} 
          isLoading={isScanning} 
        />

        {results.length > 0 && (
          <div className="space-y-10 animate-in fade-in slide-in-from-bottom-4 duration-1000">
            <SkyChart data={results} />
            <AnalysisPanel analysis={analysis} isAnalyzing={isAnalyzing} />
            
            <section className="bg-slate-950/80 border border-slate-800 rounded-2xl overflow-hidden shadow-2xl backdrop-blur-lg">
              <div className="p-5 border-b border-slate-800 bg-slate-900/40 flex justify-between items-center">
                <h3 className="font-bold text-xs tracking-widest text-slate-400 uppercase flex items-center gap-2">
    <div className="w-2 h-2 rounded-full bg-cyan-500 animate-ping"></div>
    Real-time Coordinate Inversion Log (Sorted by Depth)
  </h3>
  
  <p className="mt-1 text-[10px] italic text-slate-500 tracking-tight leading-tight max-w-md">
    NOTICE: This UI is a research prototype. Values may exhibit minor numerical 
    drift (±10%) from the stable python core due to floating-point variance.
  </p>
              </div>
              <div className="overflow-x-auto">
                <table className="w-full text-left text-xs">
                  <thead className="bg-slate-900/20 text-slate-500 uppercase tracking-tighter">
                    <tr>
                      <th className="px-6 py-4 font-semibold">XYZ Depth</th>
                      <th className="px-6 py-4 font-semibold text-cyan-400">HMS / DMS Coordinate</th>
                      <th className="px-6 py-4 font-semibold text-pink-500/80">Metric Pull</th>
                      <th className="px-6 py-4 font-semibold text-indigo-400">Sync</th>
                      <th className="px-6 py-4 font-semibold">Floor</th>
                      <th className="px-6 py-4 font-semibold">SIMBAD RESOLUTION</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-800/50">
                    {results.map((h) => (
                      <tr key={h.xyz} className="hover:bg-cyan-500/[0.03] transition-colors group">
                        <td className="px-6 py-4 text-slate-500 group-hover:text-slate-300 font-mono">
                          {h.xyz >= 1000000 ? (h.xyz / 1000000).toFixed(2) + 'M' : h.xyz.toLocaleString()}
                        </td>
                        <td className="px-6 py-4 font-medium text-cyan-100 font-mono tracking-tight">{h.hms_dms}</td>
                        <td className="px-6 py-4 font-mono font-bold text-pink-400">{h.pull.toFixed(4)}</td>
                        <td className="px-6 py-4 text-indigo-300 font-mono">{h.sync.toFixed(4)}</td>
                        <td className="px-6 py-4 text-slate-500">{h.floor.toFixed(2)}</td>
                        <td className="px-6 py-4">
                          <span className={`px-2 py-0.5 rounded text-[9px] uppercase font-bold tracking-widest border transition-all ${
                            h.simbad_id.includes('ANOMALY') 
                              ? 'bg-red-500/20 text-red-400 border-red-500/50 shadow-[0_0_10px_rgba(239,68,68,0.2)]'
                              : h.simbad_id.includes('VOID') || h.simbad_id.includes('DEEP')
                                ? 'bg-indigo-500/10 text-indigo-400 border-indigo-500/30' 
                                : 'bg-slate-800 text-slate-400 border-slate-700'
                          }`}>
                            {h.simbad_id}
                          </span>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </section>
          </div>
        )}
      </main>

      <footer className="mt-20 py-12 border-t border-slate-800 text-center opacity-40">
        <p className="text-slate-500 text-[9px] uppercase tracking-[0.4em]">
          &copy; 2025 DARK SECTOR ANALYTICS &bull; BOOTES VOID SURVEY MODE
        </p>
      </footer>
    </div>
  );
};

export default App;
