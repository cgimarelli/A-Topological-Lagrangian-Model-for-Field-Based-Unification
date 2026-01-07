
import React from 'react';
import { ScatterChart, Scatter, XAxis, YAxis, ZAxis, Tooltip, CartesianGrid, ResponsiveContainer, Cell } from 'recharts';
import { ScanResult } from '../types';

interface Props {
  data: ScanResult[];
}

const SkyChart: React.FC<Props> = ({ data }) => {
  if (data.length === 0) return null;

  const CustomTooltip = ({ active, payload }: any) => {
    if (active && payload && payload.length) {
      const hit = payload[0].payload as ScanResult;
      return (
        <div className="bg-slate-900 border border-cyan-500 p-3 rounded shadow-xl text-xs space-y-1">
          <p className="font-bold text-cyan-400 underline">{hit.simbad_id}</p>
          <p><span className="text-slate-400">Coord:</span> {hit.hms_dms}</p>
          <p><span className="text-slate-400">Pull:</span> {hit.pull.toFixed(4)}</p>
          <p><span className="text-slate-400">Sync:</span> {hit.sync.toFixed(4)}</p>
          <p><span className="text-slate-400">XYZ:</span> {hit.xyz.toFixed(2)}</p>
        </div>
      );
    }
    return null;
  };

  return (
    <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-8">
      <div className="bg-slate-900/50 border border-slate-700 p-6 rounded-xl backdrop-blur-sm h-[400px]">
        <h3 className="text-sm font-bold text-slate-400 mb-4 flex items-center gap-2 uppercase">
          <svg xmlns="http://www.w3.org/2000/svg" className="h-4 w-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M3.055 11H5a2 2 0 012 2v1a2 2 0 002 2 2 2 0 012 2v2.945M8 3.935V5.5A2.5 2.5 0 0010.5 8h.5a2 2 0 012 2 2 2 0 002 2h1.5a3.375 3.375 0 013.375 3.375V15" />
          </svg>
          Metric Tension Map (RA vs DEC)
        </h3>
        <ResponsiveContainer width="100%" height="90%">
          <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
            <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
            <XAxis 
              type="number" 
              dataKey="ra" 
              name="RA" 
              unit="h" 
              domain={['dataMin - 0.5', 'dataMax + 0.5']} 
              stroke="#64748b"
              tick={{ fontSize: 10 }}
            />
            <YAxis 
              type="number" 
              dataKey="dec" 
              name="Dec" 
              unit="°" 
              domain={['dataMin - 5', 'dataMax + 5']} 
              stroke="#64748b"
              tick={{ fontSize: 10 }}
            />
            <ZAxis type="number" dataKey="pull" range={[40, 600]} />
            <Tooltip content={<CustomTooltip />} />
            <Scatter name="Gravity Hits" data={data}>
              {data.map((entry, index) => (
                <Cell 
                  key={`cell-${index}`} 
                  fill={entry.pull > 3.0 ? '#f43f5e' : entry.pull > 1.0 ? '#f59e0b' : '#06b6d4'} 
                />
              ))}
            </Scatter>
          </ScatterChart>
        </ResponsiveContainer>
      </div>

      <div className="bg-slate-900/50 border border-slate-700 p-6 rounded-xl backdrop-blur-sm h-[400px]">
        <h3 className="text-sm font-bold text-slate-400 mb-4 flex items-center gap-2 uppercase">
          <svg xmlns="http://www.w3.org/2000/svg" className="h-4 w-4" fill="none" viewBox="0 0 24 24" stroke="currentColor">
            <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 7h8m0 0v8m0-8l-8 8-4-4-6 6" />
          </svg>
          Sync Stability vs. XYZ Distance
        </h3>
        <ResponsiveContainer width="100%" height="90%">
          <ScatterChart margin={{ top: 20, right: 20, bottom: 20, left: 20 }}>
            <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
            <XAxis 
              type="number" 
              dataKey="xyz" 
              name="XYZ Index" 
              stroke="#64748b"
              tick={{ fontSize: 10 }}
              domain={['dataMin', 'dataMax']}
            />
            <YAxis 
              type="number" 
              dataKey="sync" 
              name="Sync" 
              domain={[0.85, 1.0]} 
              stroke="#64748b"
              tick={{ fontSize: 10 }}
            />
            <Tooltip content={<CustomTooltip />} />
            <Scatter line lineJointType="monotoneX" shape="circle" data={data} fill="#8b5cf6" />
          </ScatterChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
};

export default SkyChart;
