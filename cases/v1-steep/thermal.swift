// Prints macOS's thermal state: nominal, fair, serious or critical. No admin rights needed.
import Foundation

print(["nominal", "fair", "serious", "critical"][ProcessInfo.processInfo.thermalState.rawValue])
