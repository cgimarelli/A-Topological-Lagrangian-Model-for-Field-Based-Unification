
import React from 'react';
import { ScannerParams } from '../types';

interface Props {
  params: ScannerParams;
  onChange: (params: ScannerParams) => void;
  onScan: () => void;
  isLoading: boolean;
}

const ScannerForm: React.FC<Props> = ({ params, onChange, onScan, isLoading }) => {
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    onChange({ ...params, [name]: parseFloat(value) || 0 });
  };

  return (
    <div className="bg-slate-900/50 border border-slate-700 p-6 rounded-xl backdrop-blur-sm shadow-xl">
      <h2 className="text-xl font-bold mb-6 text-cyan-400 flex items-center gap-2">
        <svg xmlns="http://www.w3.org/2000/svg" className="h-6 w-6" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 6V4m0 2a2 2 0 100 4m0-4a2 2 0 110 4m-6 8a2 2 0 100-4m0 4a2 2 0 110-4m0 4v2m0-6V4m6 6v10m6-2a2 2 0 100-4m0 4a2 2 0 110-4m0 4v2m0-6V4" />
        </svg>
        Scanner Configuration
      </h2>
      
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <div className="space-y-4">
          <div>
            <label className="block text-xs text-slate-400 mb-1">BASE RA (HOURS)</label>
            <input
              type="number"
              name="base_ra"
              step="0.1"
              value={params.base_ra}
              onChange={handleChange}
              className="w-full bg-slate-950 border border-slate-700 rounded p-2 text-cyan-50 focus:border-cyan-500 outline-none transition-colors"
            />
          </div>
          <div>
            <label className="block text-xs text-slate-400 mb-1">BASE DEC (DEGREES)</label>
            <input
              type="number"
              name="base_dec"
              step="0.1"
              value={params.base_dec}
              onChange={handleChange}
              className="w-full bg-slate-950 border border-slate-700 rounded p-2 text-cyan-50 focus:border-cyan-500 outline-none transition-colors"
            />
          </div>
        </div>
        
        <div className="space-y-4">
          <div>
            <label className="block text-xs text-slate-400 mb-1">START XYZ</label>
            <input
              type="number"
              name="start_xyz"
              value={params.start_xyz}
              onChange={handleChange}
              className="w-full bg-slate-950 border border-slate-700 rounded p-2 text-cyan-50 focus:border-cyan-500 outline-none transition-colors"
            />
          </div>
          <div>
            <label className="block text-xs text-slate-400 mb-1">SEARCH RADIUS</label>
            <input
              type="number"
              name="search_radius"
              value={params.search_radius}
              onChange={handleChange}
              className="w-full bg-slate-950 border border-slate-700 rounded p-2 text-cyan-50 focus:border-cyan-500 outline-none transition-colors"
            />
          </div>
        </div>
      </div>
      
      <button
        onClick={onScan}
        disabled={isLoading}
        className="w-full mt-8 bg-cyan-600 hover:bg-cyan-500 disabled:bg-slate-700 text-white font-bold py-3 px-6 rounded-lg shadow-lg shadow-cyan-900/20 transition-all active:scale-95 flex items-center justify-center gap-2"
      >
        {isLoading ? (
          <span className="animate-pulse">SCANNING DIMENSIONS...</span>
        ) : (
          <>
            <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
            </svg>
            INITIATE SCAN
          </>
        )}
      </button>
    </div>
  );
};

export default ScannerForm;
