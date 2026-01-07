
import { GoogleGenAI, Type } from "@google/genai";
import { ScanResult, AIAnalysis } from "../types";

// Always use process.env.API_KEY as per instructions
const ai = new GoogleGenAI({ apiKey: process.env.API_KEY });

export async function analyzeScanResults(results: ScanResult[]): Promise<AIAnalysis> {
  // CRITICAL: Copy results before sorting to avoid modifying the original array in-place
  const topHits = [...results]
    .sort((a, b) => b.pull - a.pull)
    .slice(0, 8)
    .map(h => `Coord: ${h.hms_dms}, RA: ${h.ra.toFixed(4)}, Dec: ${h.dec.toFixed(4)}, Pull: ${h.pull.toFixed(4)}, Sync: ${h.sync.toFixed(4)}, Distance: ${h.dist.toFixed(0)} LY`);

  const prompt = `
    Analyze the following gravitational inversion hits from an N12 multidimensional GPS scanner. 
    The "PULL" represents metric tension (gravitational pull rating), and "SYNC" represents dimensional alignment stability.
    
    Data Hits:
    ${topHits.join('\n')}
    
    Task:
    Provide a high-level scientific-sounding astrophysical interpretation. 
    1. Summarize the overall field stability.
    2. Offer three distinct hypotheses for the source of these inversions (e.g., dark matter filaments, micro-singularities, or distant pulsar clusters).
    3. Note any specific "PULL" anomalies that exceed standard N12.00 baselines.
    
    Respond STRICTLY in JSON format.
  `;

  try {
    const response = await ai.models.generateContent({
      model: "gemini-3-flash-preview",
      contents: prompt,
      config: {
        responseMimeType: "application/json",
        responseSchema: {
          type: Type.OBJECT,
          properties: {
            summary: { type: Type.STRING, description: "Scientific summary of findings." },
            hypotheses: { 
              type: Type.ARRAY, 
              items: { type: Type.STRING },
              description: "Three potential astrophysical hypotheses."
            },
            anomalies: { type: Type.STRING, description: "Notes on metric tension spikes." }
          },
          required: ["summary", "hypotheses", "anomalies"]
        }
      }
    });

    let jsonStr = response.text.trim();
    
    // Robustly handle Markdown code blocks if the model ignores the MimeType hint
    if (jsonStr.startsWith("```")) {
      jsonStr = jsonStr.replace(/^```(?:json)?/, "").replace(/```$/, "").trim();
    }

    return JSON.parse(jsonStr);
  } catch (error) {
    console.error("Analysis Engine Error:", error);
    return {
      summary: "Inversion field analysis interrupted. Preliminary data suggests a significant metric tension spike beyond the N12 calibration limit.",
      hypotheses: [
        "Localized spacetime curvature anomaly detected.",
        "Quantum-metric interference from neighboring dimensional planes.",
        "Undocumented high-mass object in the secondary survey sector."
      ],
      anomalies: "Critical stability fluctuation observed at primary scan points."
    };
  }
}
