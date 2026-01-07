
import React from 'react';
import { AIAnalysis } from '../types';

interface Props {
  analysis: AIAnalysis | null;
  isAnalyzing: boolean;
}

const AnalysisPanel: React.FC<Props> = ({ analysis, isAnalyzing }) => {
  if (isAnalyzing) {
    return (
      <div className="bg-cyan-950/20 border border-cyan-800 p-8 rounded-xl animate-pulse mb-8">
        <div className="h-4 bg-cyan-900/50 rounded w-1/4 mb-4"></div>
        <div className="h-20 bg-cyan-900/50 rounded mb-4"></div>
        <div className="h-4 bg-cyan-900/50 rounded w-1/3"></div>
      </div>
    );
  }

  if (!analysis) return null;

  return (
    <div className="bg-slate-900/50 border border-indigo-500/30 p-8 rounded-xl backdrop-blur-md mb-8 shadow-2xl relative overflow-hidden">
      <div className="absolute top-0 right-0 p-4 opacity-10 pointer-events-none">
        <svg xmlns="http://www.w3.org/2000/svg" className="h-32 w-32" fill="none" viewBox="0 0 24 24" stroke="currentColor">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.618.309a6 6 0 01-3.86.517l-3.158-.632a2 2 0 01-1.223-1.223L2.79 5.29a2 2 0 011.223-1.223l7-1.4a2 2 0 011.223 1.223l.542 2.711a2 2 0 001.022.547l2.387.477a6 6 0 003.86-.517l.618-.309a6 6 0 013.86-.517l3.158.632a2 2 0 011.223 1.223l.542 2.711a2 2 0 001.022.547l2.387.477a6 6 0 003.86-.517l.618-.309a6 6 0 013.86-.517l3.158.632a2 2 0 011.223 1.223l.542 2.711" />
        </svg>
      </div>

      <div className="relative z-10">
        <h2 className="text-xl font-bold text-indigo-400 mb-4 flex items-center gap-2">
          <span className="bg-indigo-900/50 p-2 rounded-lg text-indigo-300">
            AI
          </span>
          Astrophysical Interpretation
        </h2>
        
        <p className="text-slate-200 leading-relaxed mb-6 text-sm italic border-l-4 border-indigo-500 pl-4 bg-indigo-500/5 py-2">
          {analysis.summary}
        </p>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
          <div>
            <h3 className="text-xs font-bold text-slate-400 mb-3 uppercase tracking-widest">Hypotheses</h3>
            <ul className="space-y-2">
              {analysis.hypotheses.map((h, i) => (
                <li key={i} className="text-sm text-slate-300 flex items-start gap-2">
                  <span className="text-indigo-500 mt-1">•</span>
                  {h}
                </li>
              ))}
            </ul>
          </div>
          <div>
            <h3 className="text-xs font-bold text-slate-400 mb-3 uppercase tracking-widest">Detected Anomalies</h3>
            <p className="text-sm text-pink-400/80 bg-pink-900/10 p-3 rounded border border-pink-900/30">
              {analysis.anomalies}
            </p>
          </div>
        </div>
      </div>
    </div>
  );
};

export default AnalysisPanel;
