#ifndef AI_H
#define AI_H

#include <iostream>
#include <vector>
#include "Board.h"

class AI {
public:
	AI() = delete;
	static bool gameWon(Board *board, std::vector<int> &bitBoard);
	static bool gameWon(Board *board, bool isBlack);
	static Ply computeScore(Board *board, int8_t maxSearchDepth, double maxSearchTime, bool isBlack, bool printAIInformation);
private:
	static int16_t computeScoreUsingAlphaBeta(Board *board, int8_t maxLevels, int8_t levelsToBottom, int16_t alpha, int16_t beta, bool black, Ply *bestPly);
	static void alphaBetaInnerLoop(Board *board, Ply *bestPly, Position position, int16_t *alpha, int16_t *beta, int8_t levelsToBottom, int8_t maxLevels, bool isBlack, Ply *bestNextPly, bool *hasFoundABestNextPly);
	static int16_t evaluateBoard(Board *board, bool isBlack);
	static int numberOfFoursOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfFoursHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfThreesOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfThreesHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfTwosOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
	static int numberOfTwosHalfOpen(Board *board, std::vector<int> &bitBoard, std::vector<int> &bitBoardEmpty);
};

#endif
